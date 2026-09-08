process MERGE_QUANT_ENCYCLOPEDIA_DIA {
    label 'python'
    tag   "merge_enc_dia"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/base_quant/enc_dia" },
               mode: params.publish_mode

    input:
    path peptides_tsv, stageAs: 'peptides/*'
    path elib
    path enc_targets, stageAs: 'enc/*'
    path diathem_tsv
    path sample_map
    val  heavy_label
    val  primary

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

    PEP_TSV    = "${peptides_tsv}"
    ELIB       = "${elib}"
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

    _flank_re = re.compile(r'^[A-Z\\-]\\.|\\.[A-Z\\-]\$')
    def strip_mods(s):
        s = re.sub(r'\\.\\([^)]+\\)|\\([^)]+\\)|\\[[^\\]]*\\]', '', str(s))
        return _flank_re.sub('', s).strip('.')

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
    channel_of = lambda p: 'heavy' if is_heavy(p) else 'light'

    _psm_chg_re = re.compile(r'\\+(\\d+)\\s*\$')
    def charge_from_psmid(v):
        m = _psm_chg_re.search(str(v))
        return int(m.group(1)) if m else 0

    def dedupe_accessions(v):
        s = str(v or '').strip()
        if not s: return ''
        parts = [p.strip() for p in re.split(r'[;,]', s) if p.strip()]
        return ';'.join(sorted(set(parts)))

    sm = pd.read_csv(SAMPLE_MAP, sep='\\t', header=None, names=['mzml_stem','sample_id'])
    stem2sid = dict(zip(sm['mzml_stem'], sm['sample_id']))
    all_sample_ids = sm['sample_id'].tolist()

    def sid_from_path(p, sfx=''):
        stem = re.sub(re.escape(sfx) + r'\$', '', os.path.basename(p))
        for _, r in sm.iterrows():
            if stem.startswith(r['mzml_stem']) or r['mzml_stem'].startswith(stem):
                return r['sample_id']
        return stem

    quant_rows = []
    _elib_ok = os.path.exists(ELIB) and os.path.getsize(ELIB)
    if not _elib_ok and os.path.basename(PEP_TSV) != 'NO_FILE' and os.path.getsize(PEP_TSV):
        w = pd.read_csv(PEP_TSV, sep='\\t')
        meta_cols = [c for c in ('Peptide','Protein','numFragments') if c in w.columns]
        sample_cols = [c for c in w.columns if c not in meta_cols]
        for c in sample_cols:
            sid = sid_from_path(c)
            sub = w[meta_cols + [c]].rename(columns={c: 'abundance_encyclopedia'})
            sub['sample_id']    = sid
            sub['peptidoform']  = sub.get('Peptide', '')
            sub['stripped_seq'] = sub.get('Peptide', '').astype(str).apply(strip_mods)
            sub['protein']      = sub.get('Protein', '').astype(str).apply(dedupe_accessions)
            sub['charge']       = 0
            sub['channel']      = sub['peptidoform'].apply(channel_of)
            quant_rows.append(sub[['sample_id','peptidoform','stripped_seq','charge',
                                    'channel','protein','abundance_encyclopedia']])
    if not quant_rows and _elib_ok:
        import sqlite3
        con = sqlite3.connect(ELIB)
        try:
            pq = pd.read_sql_query("SELECT PeptideModSeq, PrecursorCharge, SourceFile, "
                                   "TotalIntensity FROM peptidequants", con)
            try:
                p2p = pd.read_sql_query("SELECT PeptideSeq, ProteinAccession FROM "
                                        "peptidetoprotein WHERE isDecoy = 0", con)
                pep2prot = (p2p.drop_duplicates().groupby('PeptideSeq')['ProteinAccession']
                              .apply(lambda x: ';'.join(sorted(set(x)))).to_dict())
            except Exception:
                pep2prot = {}
        finally:
            con.close()
        if not pq.empty:
            pq['peptidoform']  = pq['PeptideModSeq'].astype(str)
            pq['stripped_seq'] = pq['peptidoform'].apply(strip_mods)
            pq['channel']      = pq['peptidoform'].apply(channel_of)
            pq['sample_id']    = pq['SourceFile'].astype(str).apply(sid_from_path)
            pq = pq.rename(columns={'TotalIntensity': 'abundance_encyclopedia'})
            pq['charge'] = pd.to_numeric(pq['PrecursorCharge'], errors='coerce').fillna(0).astype(int)
            enc_fb = (pq.groupby(['sample_id','stripped_seq','channel','charge'], as_index=False)
                        .agg(abundance_encyclopedia=('abundance_encyclopedia','max'),
                             peptidoform=('peptidoform','first')))
            enc_fb['protein'] = (enc_fb['stripped_seq'].map(pep2prot).fillna('')
                                 .apply(dedupe_accessions))
            quant_rows.append(enc_fb[['sample_id','peptidoform','stripped_seq','charge',
                                      'channel','protein','abundance_encyclopedia']])
            print(f"abundance read from {ELIB} peptidequants: {len(enc_fb)} rows "
                  f"(charge-resolved)", file=sys.stderr)

    enc = pd.concat(quant_rows, ignore_index=True) if quant_rows else pd.DataFrame(
        columns=['sample_id','peptidoform','stripped_seq','charge','channel',
                 'protein','abundance_encyclopedia'])

    id_rows = []
    for pth in sorted(glob.glob('enc/*.encyclopedia.txt')):
        sid = sid_from_path(pth, '.encyclopedia.txt')
        try:
            t = pd.read_csv(pth, sep='\\t')
        except Exception as e:
            print(f"skip {pth}: {e}", file=sys.stderr); continue
        if t.empty: continue
        pep  = next((c for c in t.columns
                     if c.lower() in ('peptidemodseq','modifiedsequence','peptide','peptidoform')), None)
        chg  = next((c for c in t.columns
                     if c.lower() in ('precursorcharge','charge')), None)
        q    = next((c for c in t.columns if c.lower() in ('qvalue','q-value','q_value')), None)
        perr = next((c for c in t.columns if c.lower() in ('posterior_error_prob','posterior_error_probability','pep')), None)
        if pep is None: continue
        id_rows.append(pd.DataFrame({
            'sample_id':    sid,
            'peptidoform':  t[pep].astype(str),
            'stripped_seq': t[pep].astype(str).apply(strip_mods),
            'charge':       (t[chg].astype(int) if chg
                             else t[t.columns[0]].apply(charge_from_psmid)),
            'id_qvalue':    t[q] if q else pd.NA,
            'id_pep':       t[perr] if perr else pd.NA,
        }).assign(channel=lambda d: d['peptidoform'].apply(channel_of)))
    ids = pd.concat(id_rows, ignore_index=True) if id_rows else pd.DataFrame(
        columns=['sample_id','peptidoform','stripped_seq','charge','id_qvalue',
                 'id_pep','channel'])

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

    enc_and_ids = enc.merge(ids, on=['sample_id','stripped_seq','charge','channel'],
                             how='outer', suffixes=('_enc','_id'))
    if 'peptidoform_id' in enc_and_ids.columns:
        enc_and_ids['peptidoform'] = enc_and_ids['peptidoform_enc'].fillna(enc_and_ids['peptidoform_id'])
        enc_and_ids = enc_and_ids.drop(columns=['peptidoform_enc','peptidoform_id'])
    if 'channel' not in enc_and_ids.columns:
        enc_and_ids['channel'] = enc_and_ids['peptidoform'].apply(channel_of)

    base = enc_and_ids.merge(dia,
        on=['sample_id','stripped_seq','charge','channel'], how='outer',
        suffixes=('','_d'))
    if 'peptidoform_d' in base.columns:
        base['peptidoform'] = base['peptidoform'].fillna(base['peptidoform_d'])
        base = base.drop(columns=['peptidoform_d'])
    if 'peptide' not in base.columns:
        base['peptide'] = base['peptidoform']

    if pep2prot:
        _blank = base['protein'].isna() | base['protein'].astype(str).str.strip().eq('')
        base.loc[_blank, 'protein'] = (base.loc[_blank, 'stripped_seq']
                                         .map(pep2prot).fillna(''))
    for c in ['abundance_openswath','abundance_diathem','abundance_encyclopedia',
              'abundance_tric','consistency','n_effective_transitions','rt_apex_seconds']:
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
            'id_qvalue','id_pep','calibrated']
    cols = [c for c in cols if c in base.columns]
    out = base[cols].sort_values(['sample_id','peptide','charge','channel'])
    out.to_csv(OUT_LONG, sep='\\t', index=False)
    for sid in all_sample_ids:
        grp = out[out['sample_id'] == sid]
        grp.to_csv(f"{OUT_DIR}/{sid}.base_quant.tsv", sep='\\t', index=False)

    print(f"enc rows: {len(enc)}", file=sys.stderr)
    print(f"id rows:  {len(ids)}", file=sys.stderr)
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
