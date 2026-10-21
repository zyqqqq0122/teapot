process MERGE_QUANT_OPENSWATH_PRM {
    label 'python'
    tag   "merge_osw_prm"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/base_quant/osw_prm" },
               mode: params.publish_mode

    input:
    path context_peptide_tables, stageAs: 'ids/*'
    path osw_intensities,    stageAs: 'osw_int/*'
    path diathem_tsv
    path sample_map
    val  heavy_label
    val  primary

    output:
    path "base_quant.long.tsv",         emit: long
    path "base_quant_by_sample/*.tsv",  emit: per_sample
    path "merge_quant.log",             emit: log
    path "versions.yml",                emit: versions

    script:
    def _rows_uni = file("${projectDir}/assets/unimod_delta_masses.tsv").readLines()
    def _delta_by_acc = [:]
    _rows_uni.drop(1).each { l ->
        def c = l.split('\t')
        if (c.size() >= 3) _delta_by_acc[c[0]] = c[2] as Double
    }
    def _aug = (heavy_label ?: []).collect { r ->
        def d = r.containsKey('delta') ? (r.delta as Double) : _delta_by_acc[r.unimod]
        [ site: r.site, unimod: r.unimod, delta: d ]
    }
    def label_json = groovy.json.JsonOutput.toJson(_aug)
    def mod_tol    = params.mod_mass_tolerance
    """
    set -o pipefail

    python3 - <<'PY' 2>&1 | tee merge_quant.log
    import pandas as pd, glob, os, re, json, sys

    DIA_TSV    = "${diathem_tsv}"
    SAMPLE_MAP = "${sample_map}"
    LABEL      = json.loads(r'''${label_json}''')
    MTOL       = float("${mod_tol}")
    PRIMARY    = "${primary}"
    OUT_LONG   = "base_quant.long.tsv"
    OUT_DIR    = "base_quant_by_sample"
    os.makedirs(OUT_DIR, exist_ok=True)

    _accs   = [r["unimod"] for r in LABEL if r.get("unimod")]
    _deltas = [float(r["delta"]) for r in LABEL if r.get("delta") is not None]
    _acc_re   = re.compile("|".join(re.escape(a) for a in _accs)) if _accs else None
    _delta_re = re.compile(r'[\\[\\(]([+-]?\\d+(?:\\.\\d+)?)[\\]\\)]')
    HAS_HEAVY = bool(LABEL)

    def strip_mods(s):
        return re.sub(r'\\.\\([^)]+\\)|\\([^)]+\\)|\\[[^\\]]*\\]', '', str(s)).strip('.')

    def is_heavy(peptidoform):
        if not HAS_HEAVY:
            return False
        s = str(peptidoform)
        if _acc_re and _acc_re.search(s):
            return True
        if _deltas:
            for m in _delta_re.finditer(s):
                try:
                    v = float(m.group(1))
                except ValueError:
                    continue
                for d in _deltas:
                    if abs(v - d) <= MTOL:
                        return True
        return False

    def channel_of(peptidoform):
        return 'heavy' if is_heavy(peptidoform) else 'light'

    sm = pd.read_csv(SAMPLE_MAP, sep='\\t', header=None,
                     names=['mzml_stem', 'sample_id'])
    all_sample_ids = sm['sample_id'].tolist()

    osw_rows = []
    for pth in sorted(glob.glob('osw_int/*.osw_intensity.tsv')):
        sid = re.sub(r'\\.osw_intensity\\.tsv\$', '', os.path.basename(pth))
        try:
            t = pd.read_csv(pth, sep='\\t')
        except Exception as e:
            print(f"skip {pth}: {e}", file=sys.stderr); continue
        if t.empty:
            continue
        t = t[t.get('decoy', 0).astype(int) == 0]
        t['sample_id']   = sid
        t['channel']     = t['peptidoform'].apply(channel_of)
        osw_rows.append(t[['sample_id','stripped_seq','peptidoform','charge',
                           'channel','rt_apex_seconds','abundance_openswath'] +
                          [c for c in ['width_openswath_seconds'] if c in t.columns]])
    osw = pd.concat(osw_rows, ignore_index=True) if osw_rows else pd.DataFrame(
        columns=['sample_id','stripped_seq','peptidoform','charge','channel',
                 'rt_apex_seconds','abundance_openswath'])

    id_rows = []
    for pep_path in sorted(glob.glob('ids/*.peptide.reference.txt')):
        sid = re.sub(r'\\.peptide\\.reference\\.txt\$', '', os.path.basename(pep_path))
        try:
            t = pd.read_csv(pep_path, sep='\\t')
        except Exception as e:
            print(f"skip {pep_path}: {e}", file=sys.stderr); continue
        if t.empty:
            continue
        pep_col  = next((c for c in t.columns
                         if c.lower() in ('peptide', 'sequence', 'peptidesequence')), t.columns[-2])
        prot_col = 'Proteins' if 'Proteins' in t.columns else t.columns[-1]
        q_col    = next((c for c in t.columns if c.lower() in ('q-value','q_value','qvalue')), None)
        rank_col = next((c for c in t.columns if c.lower() == 'peak_group_rank'), None)
        perr_col = next((c for c in t.columns if c.lower() in ('posterior_error_prob','posterior_error_probability','pep')), None)
        id_rows.append(pd.DataFrame({
            'sample_id':    sid,
            'peptide':      t[pep_col].astype(str),
            'stripped_seq': t[pep_col].astype(str).apply(strip_mods),
            'channel':      t[pep_col].astype(str).apply(channel_of),
            'protein':      t[prot_col].astype(str),
            'id_qvalue':    t[q_col] if q_col else pd.NA,
            'id_pep':       t[perr_col] if perr_col else pd.NA,
            'peak_group_rank': t[rank_col] if rank_col else pd.NA,
        }))
    ids = pd.concat(id_rows, ignore_index=True) if id_rows else pd.DataFrame(
        columns=['sample_id','peptide','stripped_seq','channel','protein','id_qvalue',
                 'id_pep','peak_group_rank'])

    if os.path.basename(DIA_TSV) == 'NO_FILE' or not os.path.getsize(DIA_TSV):
        dia = pd.DataFrame(columns=['sample_id','peptidoform','charge','channel',
                                    'stripped_seq','abundance_diathem',
                                    'consistency','n_effective_transitions',
                                    'rt_diathem_seconds'])
    else:
        dia = pd.read_csv(DIA_TSV, sep='\\t')
        if 'is_decoy' in dia.columns:
            dia = dia[pd.to_numeric(dia['is_decoy'], errors='coerce').fillna(0) == 0]
        dia = dia.rename(columns={'run_id': 'sample_id',
                                   'abundance': 'abundance_diathem',
                                   'rt_apex_seconds': 'rt_diathem_seconds'})
        dia['stripped_seq'] = dia['peptidoform'].apply(strip_mods)
        dia['channel']      = dia['peptidoform'].apply(channel_of)
        if 'rt_diathem_seconds' not in dia.columns: dia['rt_diathem_seconds'] = pd.NA
        dia = dia[['sample_id','peptidoform','charge','channel','stripped_seq',
                   'abundance_diathem','consistency','n_effective_transitions',
                   'rt_diathem_seconds']]

    key = ['sample_id','stripped_seq','charge','channel']
    base = osw.merge(dia, on=key, how='outer', suffixes=('_osw','_dia'))
    base['peptidoform'] = base['peptidoform_osw'].fillna(base['peptidoform_dia'])
    base = base.drop(columns=[c for c in ['peptidoform_osw','peptidoform_dia'] if c in base.columns])

    base = base.merge(ids, on=['sample_id','stripped_seq','channel'], how='outer')
    base['peptide'] = base['peptide'].fillna(base['peptidoform'])

    for c in ['abundance_openswath','abundance_diathem','abundance_encyclopedia',
              'abundance_tric','consistency','n_effective_transitions',
              'rt_apex_seconds','rt_diathem_seconds','width_openswath_seconds']:
        if c not in base.columns: base[c] = pd.NA
    base['channel'] = base['channel'].fillna('light')

    METHODS = ['encyclopedia','openswath','diathem','tric']
    _pcol = f'abundance_{PRIMARY}'
    base['abundance_primary'] = base[_pcol] if _pcol in base.columns else pd.NA
    base['abundance_primary_source'] = base['abundance_primary'].map(
        lambda v: PRIMARY if pd.notna(v) else None)
    _n_pri = int(base['abundance_primary'].notna().sum())
    print(f"abundance_primary: {PRIMARY} only, no cross-tool fallback "
          f"({_n_pri} of {len(base)} rows have a value)", file=sys.stderr)

    if HAS_HEAVY:
        by_key = base.groupby(['sample_id','stripped_seq','charge'])['channel'] \\
                     .apply(lambda ch: 'light' in set(ch) and 'heavy' in set(ch))
        base['calibrated'] = base.set_index(['sample_id','stripped_seq','charge']) \\
                                 .index.map(by_key).fillna(False).astype(bool)
    else:
        base['calibrated'] = False

    cols = ['sample_id','peptide','peptidoform','charge','channel','protein',
            'abundance_openswath','abundance_encyclopedia','abundance_diathem',
            'abundance_tric','abundance_primary','abundance_primary_source',
            'consistency','n_effective_transitions','rt_apex_seconds','rt_diathem_seconds',
            'width_openswath_seconds','id_qvalue','id_pep','peak_group_rank','calibrated']
    cols = [c for c in cols if c in base.columns]
    out = base[cols].sort_values(['sample_id','peptide','charge','channel'])
    out.to_csv(OUT_LONG, sep='\\t', index=False)
    for sid in all_sample_ids:
        grp = out[out['sample_id'] == sid]
        grp.to_csv(f"{OUT_DIR}/{sid}.base_quant.tsv", sep='\\t', index=False)

    print(f"osw rows:     {len(osw)}", file=sys.stderr)
    print(f"dia rows:     {len(dia)}", file=sys.stderr)
    print(f"id rows:      {len(ids)}", file=sys.stderr)
    print(f"base rows:    {len(out)}", file=sys.stderr)
    print(f"heavy rows:   {int((out['channel']=='heavy').sum())}", file=sys.stderr)
    print(f"calibratable: {int(out['calibrated'].sum())}", file=sys.stderr)
    for m in METHODS:
        c = f'abundance_{m}'
        if c in out.columns:
            print(f"  {c}: {int(out[c].notna().sum())} non-null", file=sys.stderr)
    PY

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
      primary: '${primary}'
    EOF
    """

    stub:
    """
    mkdir -p base_quant_by_sample
    touch base_quant.long.tsv merge_quant.log
    for sid in \$(cut -f2 ${sample_map}); do
        touch base_quant_by_sample/\${sid}.base_quant.tsv
    done
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
