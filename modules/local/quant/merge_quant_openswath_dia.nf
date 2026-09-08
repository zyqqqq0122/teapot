process MERGE_QUANT_OPENSWATH_DIA {
    label 'python'
    tag   "merge_osw_dia"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/base_quant/osw_dia" },
               mode: params.publish_mode

    input:
    path  pyprophet_tsvs, stageAs: 'pp/*'
    path  diathem_tsv, stageAs: 'diathem_in/*'
    path  tric_tsv,    stageAs: 'tric_in/*'
    path  sample_map
    val   heavy_label
    val   primary

    output:
    path "base_quant.long.tsv",        emit: long
    path "base_quant_by_sample/*.tsv", emit: per_sample
    path "merge_quant.log",            emit: log
    path "versions.yml",               emit: versions

    script:
    def _rows_uni = file("${projectDir}/assets/unimod_delta_masses.tsv").readLines()
    def _delta_by_acc = [:]
    _rows_uni.drop(1).each { l ->
        def c = l.split('\t')
        if (c.size() >= 3) _delta_by_acc[c[0]] = c[2] as Double
    }
    def _aug = (heavy_label ?: []).collect { r ->
        def d = r.containsKey('delta') ? (r.delta as Double) : _delta_by_acc[r.unimod]
        [ site: (r.site ?: r.residue), unimod: r.unimod, delta: d ]
    }
    def label_json = groovy.json.JsonOutput.toJson(_aug)
    def mod_tol    = params.mod_mass_tolerance
    """

    python3 - <<'PY' > merge_quant.log 2>&1
    import pandas as pd, glob, os, re, json, sys

    DIA_TSV    = "${diathem_tsv}"
    TRIC_TSV   = "${tric_tsv}"
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

    sm = pd.read_csv(SAMPLE_MAP, sep='\\t', header=None, names=['mzml_stem','sample_id'])
    all_sample_ids = sm['sample_id'].tolist()

    id_rows = []
    for pth in sorted(glob.glob('pp/*.tsv')):
        sid = re.sub(r'\\.tsv\$', '', os.path.basename(pth))
        try:
            t = pd.read_csv(pth, sep='\\t')
        except Exception as e:
            print(f"skip {pth}: {e}", file=sys.stderr); continue
        if t.empty:
            continue
        seq  = 'Sequence'        if 'Sequence'        in t.columns else 'peptide_sequence'
        full = 'FullPeptideName' if 'FullPeptideName' in t.columns else 'modified_sequence'
        chg  = 'Charge'          if 'Charge'          in t.columns else 'charge'
        inten= 'Intensity'       if 'Intensity'       in t.columns else 'intensity'
        q    = 'm_score'         if 'm_score'         in t.columns else None
        rank = 'peak_group_rank' if 'peak_group_rank' in t.columns else None
        perr = next((c for c in t.columns if c.lower() in ('posterior_error_prob','posterior_error_probability','pep')), None)
        prot = 'ProteinName'     if 'ProteinName'     in t.columns else 'protein'
        decoy= 'decoy'           if 'decoy'           in t.columns else None
        if decoy:
            t = t[t[decoy].astype(int) == 0]
        id_rows.append(pd.DataFrame({
            'sample_id':           sid,
            'peptide':             t[full].astype(str),
            'peptidoform':         t[full].astype(str),
            'stripped_seq':        t[seq].astype(str).apply(strip_mods),
            'charge':              t[chg].astype(int),
            'protein':             t[prot].astype(str),
            'abundance_openswath': t[inten] if inten in t.columns else pd.NA,
            'id_qvalue':           t[q] if q else pd.NA,
            'peak_group_rank':     t[rank] if rank else pd.NA,
            'id_pep':              t[perr] if perr else pd.NA,
        }))
    ids = pd.concat(id_rows, ignore_index=True) if id_rows else pd.DataFrame(
        columns=['sample_id','peptide','peptidoform','stripped_seq','charge',
                 'protein','abundance_openswath','id_qvalue','id_pep',
                 'peak_group_rank'])
    ids['channel'] = ids['peptidoform'].apply(channel_of)

    if os.path.basename(DIA_TSV) == 'NO_FILE' or not os.path.getsize(DIA_TSV):
        dia = pd.DataFrame(columns=['sample_id','peptidoform','charge','channel',
                                    'stripped_seq','abundance_diathem',
                                    'consistency','n_effective_transitions'])
    else:
        dia = pd.read_csv(DIA_TSV, sep='\\t')
        dia = dia.rename(columns={'run_id': 'sample_id', 'abundance': 'abundance_diathem'})
        dia['stripped_seq'] = dia['peptidoform'].apply(strip_mods)
        dia['channel']      = dia['peptidoform'].apply(channel_of)
        dia = dia[['sample_id','peptidoform','charge','channel','stripped_seq',
                   'abundance_diathem','consistency','n_effective_transitions']]

    tric_supplied = (os.path.basename(TRIC_TSV) != 'NO_FILE'
                     and os.path.getsize(TRIC_TSV) > 0)
    if not tric_supplied:
        tric = pd.DataFrame(columns=['sample_id','peptidoform','charge','channel',
                                     'abundance_tric'])
    else:
        tr = pd.read_csv(TRIC_TSV, sep='\t')
        fcol = next((c for c in ('align_origfilename','filename','align_runid','run_id')
                     if c in tr.columns), None)
        pcol = next((c for c in ('FullPeptideName','modified_sequence','Sequence')
                     if c in tr.columns), None)
        ccol = next((c for c in ('Charge','charge') if c in tr.columns), None)
        icol = next((c for c in ('Intensity','intensity') if c in tr.columns), None)
        if not all((fcol, pcol, ccol, icol)):
            sys.exit(
                "HARD ERROR: TRIC output cannot be joined -- run/file=%r peptide=%r "
                "charge=%r intensity=%r. Columns present: %s"
                % (fcol, pcol, ccol, icol, list(tr.columns)))
        dcol = next((c for c in ('decoy','Decoy') if c in tr.columns), None)
        if dcol:
            tr = tr[pd.to_numeric(tr[dcol], errors='coerce').fillna(0).astype(int) == 0]
        stem2sid = dict(zip(sm['mzml_stem'].astype(str), sm['sample_id'].astype(str)))
        def to_sid(v):
            b = os.path.basename(str(v))
            for ext in ('.mzML.gz','.mzML','.mzXML','.osw','.tsv'):
                if b.endswith(ext):
                    b = b[:-len(ext)]
                    break
            return stem2sid.get(b, b)
        tric = pd.DataFrame({
            'sample_id':      tr[fcol].apply(to_sid),
            'peptidoform':    tr[pcol].astype(str),
            'charge':         pd.to_numeric(tr[ccol], errors='coerce'),
            'abundance_tric': pd.to_numeric(tr[icol], errors='coerce'),
        }).dropna(subset=['charge'])
        tric['charge']       = tric['charge'].astype(int)
        tric['stripped_seq'] = tric['peptidoform'].apply(strip_mods)
        tric['channel']      = tric['peptidoform'].apply(channel_of)
        tric['channel'] = tric['peptidoform'].apply(channel_of)
        tric = (tric.groupby(['sample_id','peptidoform','charge','channel'],
                             dropna=False)
                    .agg(abundance_tric=('abundance_tric','sum'))
                    .reset_index())

    key = ['sample_id','stripped_seq','charge','channel']
    base = ids.merge(dia, on=key, how='outer', suffixes=('_id','_dia'))
    base['peptidoform'] = base['peptidoform_id'].fillna(base['peptidoform_dia']) \\
                            if 'peptidoform_id' in base.columns else base['peptidoform']
    base = base.drop(columns=[c for c in ['peptidoform_id','peptidoform_dia'] if c in base.columns])
    base['peptide'] = base['peptide'].fillna(base['peptidoform'])

    base = base.merge(tric, on=['sample_id','peptidoform','charge','channel'],
                      how='left')
    if tric_supplied:
        n_tric = int(base['abundance_tric'].notna().sum())
        print(f"TRIC: {len(tric)} aligned rows, {n_tric} joined onto the quant table",
              file=sys.stderr)
        if n_tric == 0:
            sys.exit(
                "HARD ERROR: TRIC ran and produced %d rows, but none joined onto the "
                "quant table on (sample_id, peptidoform, charge, channel). The likely "
                "cause is run naming: TRIC identifies runs by file, and those names "
                "must map through sample_map. Refusing to emit an all-empty "
                "abundance_tric column silently." % len(tric))

    for c in ['abundance_openswath','abundance_diathem','abundance_encyclopedia',
              'abundance_tric','consistency','n_effective_transitions',
              'rt_apex_seconds']:
        if c not in base.columns: base[c] = pd.NA

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
            'consistency','n_effective_transitions','rt_apex_seconds',
            'id_qvalue','id_pep','peak_group_rank','calibrated']
    cols = [c for c in cols if c in base.columns]
    out = base[cols].sort_values(['sample_id','peptide','charge','channel'])
    out.to_csv(OUT_LONG, sep='\\t', index=False)
    for sid in all_sample_ids:
        grp = out[out['sample_id'] == sid]
        grp.to_csv(f"{OUT_DIR}/{sid}.base_quant.tsv", sep='\\t', index=False)

    print(f"pp rows:  {len(ids)}", file=sys.stderr)
    print(f"dia rows: {len(dia)}", file=sys.stderr)
    print(f"base:     {len(out)}", file=sys.stderr)
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
