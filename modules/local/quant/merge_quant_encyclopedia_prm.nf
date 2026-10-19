process MERGE_QUANT_ENCYCLOPEDIA_PRM {
    label 'python'
    tag   "merge_enc_prm"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/base_quant/enc_prm" },
               mode: params.publish_mode

    input:
    path targets_txts,  stageAs: 'ref_tgts/*'
    path features_txts, stageAs: 'feats/*'
    path diathem_tsv
    path sample_map
    val  heavy_label
    val  primary
    path elib

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
        [ site: r.site, unimod: r.unimod, delta: d ]
    }
    def label_json = groovy.json.JsonOutput.toJson(_aug)
    def mod_tol    = params.mod_mass_tolerance
    """

    python3 - <<'PY' > merge_quant.log 2>&1
    import pandas as pd, numpy as np, glob, os, re, json, sys, functools, zlib

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
        s = _flank_re.sub('', s)
        return s.strip('.')

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

    _psm_rt_re = re.compile(r':(\\d+(?:\\.\\d*)?):')
    def rt_from_psmid(v):
        m = _psm_rt_re.search(str(v))
        return float(m.group(1)) if m else float('nan')

    def dedupe_accessions(v):
        s = str(v or '').strip()
        if not s: return ''
        parts = [p.strip() for p in re.split(r'[;,]', s) if p.strip()]
        return ';'.join(sorted(set(parts)))

    sm = pd.read_csv(SAMPLE_MAP, sep='\\t', header=None, names=['mzml_stem','sample_id'],
                     dtype=str, keep_default_na=False)
    stem2sid = dict(zip(sm['mzml_stem'], sm['sample_id']))
    all_sample_ids = sm['sample_id'].tolist()

    def sid_from_path(p, sfx=''):
        stem = re.sub(re.escape(sfx) + r'\$', '', os.path.basename(str(p)))
        for s in (stem, os.path.splitext(stem)[0]):
            if s in stem2sid:
                return stem2sid[s]
        sys.exit(f"MERGE_QUANT_ENCYCLOPEDIA_PRM: {p} is none of the runs in {SAMPLE_MAP}")

    SHARED = 'abundance_encyclopedia_shared_ions'
    PROTON, H_ATOM, WATER, NH3 = 1.00727646681290, 1.00782503207, 18.0105646837, 17.0265491
    AA_MASS = {'A': 71.037114, 'R': 156.101111, 'N': 114.042927, 'D': 115.026943,
               'C': 103.009185, 'E': 129.042593, 'Q': 128.058578, 'G': 57.021464,
               'H': 137.058912, 'I': 113.084064, 'L': 113.084064, 'K': 128.094963,
               'M': 131.040485, 'F': 147.068414, 'P': 97.052764, 'S': 87.032028,
               'T': 101.047679, 'W': 186.079313, 'Y': 163.063329, 'V': 99.068414,
               'U': 150.953636, 'O': 237.147727}
    ION_TOL = 0.005
    _res_re = re.compile(r'\\[([+-]?\\d+(?:\\.\\d*)?)\\]|([A-Z])')

    @functools.lru_cache(maxsize=None)
    def fragment_ions(modseq, charge):
        res, nterm = [], 0.0
        for mod, aa in _res_re.findall(modseq):
            if aa:
                if aa not in AA_MASS: return (), np.empty(0)
                res.append(AA_MASS[aa] + nterm); nterm = 0.0
            elif res: res[-1] += float(mod)
            else: nterm += float(mod)
        names, mzs, n = [], [], len(res)
        b, y = PROTON, WATER + PROTON
        for i in range(1, n + 1):
            b += res[i - 1]; y += res[n - i]
            for t, m in (('b', b), ('y', y), ('c', b + NH3), ('z', y - NH3),
                         ('z+1', y - NH3 + H_ATOM)):
                for z in range(1, max(1, min(charge - 1, 4)) + 1):
                    names.append((t, i, z)); mzs.append((m + (z - 1) * PROTON) / z)
        return tuple(names), np.array(mzs)

    def ion_labels(modseq, charge, mzs):
        names, theo = fragment_ions(modseq, charge)
        return [{names[k] for k in np.flatnonzero(np.abs(theo - m) <= ION_TOL)} for m in mzs]

    def shared_ion_abundance(top, keys):
        val = pd.Series(float('nan'), index=top.index)
        rows = {}
        for r in top.itertuples():
            rows.setdefault((r.sample_id, r.stripped_seq, r.charge), {})[r.channel] = r
        n_shared, same, unlabelled, no_ions = [], 0, 0, 0
        for ch in rows.values():
            l, h = ch.get('light'), ch.get('heavy')
            if l is None or h is None or not (l.abundance_encyclopedia > 0 and
                                               h.abundance_encyclopedia > 0):
                continue
            if l.QuantIonIntensityArray is None or h.QuantIonIntensityArray is None:
                no_ions += 1
                continue
            ions = [(np.frombuffer(zlib.decompress(r.QuantIonMassArray), '>f8'),
                     np.frombuffer(zlib.decompress(r.QuantIonIntensityArray), '>f4'))
                    for r in (l, h)]
            ll, hl = [ion_labels(r.PeptideModSeq, int(r.charge), mz)
                      for r, (mz, _it) in zip((l, h), ions)]
            unlabelled += sum(not s for s in ll + hl)
            li, hi = [], []
            for i, a in enumerate(ll):
                j = next((j for j, b in enumerate(hl) if j not in hi and a & b), None)
                if j is not None:
                    li.append(i); hi.append(j)
            for r, (mz, it), idx in zip((l, h), ions, (li, hi)):
                val[r.Index] = (r.abundance_encyclopedia if len(idx) == len(mz)
                                else float(it[idx].astype(float).sum()))
            n_shared.append(len(li))
            same += len(li) == len(ions[0][0]) == len(ions[1][0])
        n = pd.Series(n_shared, dtype=int)
        print(f"light/heavy pairs both quantified: {len(n)}, each form summed over the quant ions "
              f"both share (ions shared: pairs) {n.value_counts().sort_index().to_dict()}; "
              f"same ions in both forms: {same}; none shared: {int((n == 0).sum())}; quant ions "
              f"without a b/y/c/z label: {unlabelled}; pairs without ion intensities: {no_ions}",
              file=sys.stderr)
        return top[keys].assign(**{SHARED: val})

    import sqlite3
    enc_cols = ['sample_id','peptidoform','stripped_seq','charge','channel',
                'protein','abundance_encyclopedia','rt_encyclopedia_seconds',
                'width_encyclopedia_seconds'] + \\
               ([SHARED] if HAS_HEAVY else [])
    pep2prot = {}
    if os.path.exists(ELIB) and os.path.getsize(ELIB):
        con = sqlite3.connect(ELIB)
        try:
            pq = pd.read_sql_query("SELECT PeptideModSeq, PrecursorCharge, SourceFile, "
                                   "TotalIntensity, RTInSecondsCenter, RTInSecondsStart, RTInSecondsStop, "
                                   "QuantIonMassArray, "
                                   "QuantIonIntensityArray FROM peptidequants", con)
            try:
                p2p = pd.read_sql_query("SELECT PeptideSeq, ProteinAccession "
                                        "FROM peptidetoprotein WHERE isDecoy = 0", con)
                pep2prot = (p2p.drop_duplicates()
                               .groupby('PeptideSeq')['ProteinAccession']
                               .apply(lambda s: ';'.join(sorted(set(s)))).to_dict())
            except Exception as e:
                print(f"peptidetoprotein unavailable in {ELIB}: {e}", file=sys.stderr)
        finally:
            con.close()
    else:
        pq = pd.DataFrame(columns=['PeptideModSeq','PrecursorCharge','SourceFile','TotalIntensity',
                                   'RTInSecondsCenter'])
    if not pq.empty:
        pq['peptidoform']  = pq['PeptideModSeq'].astype(str)
        pq['stripped_seq'] = pq['PeptideModSeq'].astype(str).apply(strip_mods)
        pq['channel']      = pq['peptidoform'].apply(channel_of)
        pq['sample_id']    = pq['SourceFile'].astype(str).apply(sid_from_path)
        pq = pq.rename(columns={'TotalIntensity': 'abundance_encyclopedia'})
        pq['charge'] = pd.to_numeric(pq['PrecursorCharge'], errors='coerce').fillna(0).astype(int)
        keys = ['sample_id','stripped_seq','channel','charge']
        enc = (pq.groupby(keys, as_index=False)
                 .agg(abundance_encyclopedia=('abundance_encyclopedia','max'),
                      peptidoform=('peptidoform','first')))
        top = pq.loc[pq.groupby(keys)['abundance_encyclopedia'].idxmax()]
        top = top.assign(width_encyclopedia_seconds=top['RTInSecondsStop'] - top['RTInSecondsStart'])
        enc = enc.merge(top[keys + ['RTInSecondsCenter', 'width_encyclopedia_seconds']]
                        .rename(columns={'RTInSecondsCenter': 'rt_encyclopedia_seconds'}),
                        on=keys, how='left')
        if HAS_HEAVY:
            enc = enc.merge(shared_ion_abundance(top, keys), on=keys, how='left')
        enc['protein'] = enc['stripped_seq'].map(pep2prot).fillna('').apply(dedupe_accessions)
        enc = enc[enc_cols]
    else:
        enc = pd.DataFrame(columns=enc_cols)

    id_rows = []
    for pth in sorted(glob.glob('ref_tgts/*_reference.features.pep.output.txt')):
        sid = sid_from_path(pth, '_reference.features.pep.output.txt')
        try:
            t = pd.read_csv(pth, sep='\\t')
        except Exception as e:
            print(f"skip {pth}: {e}", file=sys.stderr); continue
        t = t[~t.iloc[:, 0].astype(str).str.startswith('pi_0=')]
        if t.empty: continue
        pep  = next((c for c in t.columns
                     if c.lower() in ('peptide','peptidemodseq','modifiedsequence','peptidoform')), t.columns[0])
        chg  = next((c for c in t.columns
                     if c.lower() in ('precursorcharge','charge')), None)
        q    = next((c for c in t.columns
                     if c.lower() in ('qvalue','q-value','q_value')), None)
        perr = next((c for c in t.columns
                     if c.lower() in ('posterior_error_prob','posterior_error_probability','pep')), None)
        prot = next((c for c in t.columns
                     if c.lower() in ('proteinids','proteinid','proteinaccession',
                                      'protein','proteinname','proteins',
                                      'protein_id','accession')), None)
        if prot is None:
            print(f"WARNING {pth}: no protein column among {list(t.columns)}; "
                  "absolute quant against a protein-keyed standards file "
                  "cannot resolve.", file=sys.stderr)
        psm = next((c for c in t.columns if c.lower() == 'psmid'), t.columns[0])
        d = pd.DataFrame({
            'sample_id':    sid,
            'peptidoform':  t[pep].astype(str),
            'stripped_seq': t[pep].astype(str).apply(strip_mods),
            'charge':       (t[chg].astype(int) if chg
                             else t[psm].apply(charge_from_psmid)),
            'protein':      (t[prot].astype(str).apply(dedupe_accessions)
                             if prot else ''),
            'id_qvalue':    t[q] if q else pd.NA,
            'id_pep':       t[perr] if perr else pd.NA,
            'rt_apex_seconds': t[psm].apply(rt_from_psmid),
        })
        d['channel']  = d['peptidoform'].apply(channel_of)
        d['peptide']  = d['peptidoform']
        id_rows.append(d)
    ids = pd.concat(id_rows, ignore_index=True) if id_rows else pd.DataFrame(
        columns=['sample_id','peptide','peptidoform','stripped_seq','charge','protein',
                 'id_qvalue','id_pep','rt_apex_seconds','channel'])

    if os.path.basename(DIA_TSV) == 'NO_FILE' or not os.path.getsize(DIA_TSV):
        dia = pd.DataFrame(columns=['sample_id','peptidoform','charge','channel',
                                    'stripped_seq','abundance_diathem',
                                    'consistency','n_effective_transitions',
                                    'rt_diathem_seconds'])
    else:
        dia = pd.read_csv(DIA_TSV, sep='\\t')
        if 'is_decoy' in dia.columns:
            dia = dia[pd.to_numeric(dia['is_decoy'], errors='coerce').fillna(0) == 0]
        dia = dia.rename(columns={'run_id': 'sample_id', 'abundance': 'abundance_diathem',
                                  'rt_apex_seconds': 'rt_diathem_seconds'})
        if 'rt_diathem_seconds' not in dia.columns: dia['rt_diathem_seconds'] = pd.NA
        dia['stripped_seq'] = dia['peptidoform'].apply(strip_mods)
        dia['channel']      = dia['peptidoform'].apply(channel_of)
        dia = dia[['sample_id','peptidoform','charge','channel','stripped_seq',
                   'abundance_diathem','consistency','n_effective_transitions',
                   'rt_diathem_seconds']]

    enc_and_ids = ids.merge(enc, on=['sample_id','stripped_seq','channel','charge'],
                             how='outer', suffixes=('_id','_enc'))
    for col in ('peptidoform', 'protein'):
        if f'{col}_id' in enc_and_ids.columns:
            _left = enc_and_ids[f'{col}_id'].replace('', pd.NA)
            enc_and_ids[col] = _left.fillna(enc_and_ids.get(f'{col}_enc'))
            enc_and_ids = enc_and_ids.drop(columns=[f'{col}_id', f'{col}_enc'], errors='ignore')

    base = enc_and_ids.merge(dia,
        on=['sample_id','stripped_seq','charge','channel'], how='outer',
        suffixes=('','_d'))
    if 'peptidoform_d' in base.columns:
        base['peptidoform'] = base['peptidoform'].fillna(base['peptidoform_d'])
        base = base.drop(columns=['peptidoform_d'])
    if 'peptide' not in base.columns:
        base['peptide'] = base['peptidoform']
    base['peptide'] = base['peptide'].fillna(base['peptidoform'])
    for c in ['abundance_openswath','abundance_diathem','abundance_encyclopedia',
              'abundance_tric','consistency','n_effective_transitions','rt_apex_seconds',
              'rt_encyclopedia_seconds','rt_diathem_seconds','width_encyclopedia_seconds']:
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
            'abundance_openswath','abundance_encyclopedia',SHARED,'abundance_diathem',
            'abundance_tric','abundance_primary','abundance_primary_source',
            'consistency','n_effective_transitions','rt_apex_seconds',
            'rt_encyclopedia_seconds','rt_diathem_seconds','width_encyclopedia_seconds',
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
