process QUANTIFY_HEAVY_LIGHT {
    label 'python'
    tag   "abs_quant"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/quantify_heavy_light" },
               mode: params.publish_mode

    input:
    path  base_long
    path  standard_amounts
    val   heavy_label
    val   primary
    val   min_consistency

    output:
    path "peptide.abs_quant.long.tsv",   emit: peptide_long
    path "protein.abs_quant.long.tsv",   emit: protein_long
    path "peptide.abs_quant.matrix.tsv", emit: peptide_matrix
    path "peptide.pairs.tsv",            emit: pairs
    path "peptide.pairs.confident.tsv",  emit: pairs_confident
    path "abs_quant_qc.tsv",             emit: qc
    path "quantify_heavy_light.log",        emit: log
    path "versions.yml",                 emit: versions

    script:
    def label_json = groovy.json.JsonOutput.toJson(heavy_label ?: [])
    def skip_json  = groovy.json.JsonOutput.toJson(params.complete_channels_skip_proteins ?: [])
    """

    python3 - <<'PY' > quantify_heavy_light.log 2>&1
    import pandas as pd, re, sys, json

    BASE      = "${base_long}"
    STD       = "${standard_amounts}"
    LABEL     = json.loads(r'''${label_json}''')
    PRIMARY   = "${primary}"
    MIN_CONS  = float("${min_consistency}")
    FDR       = float("${params.fdr}")
    SKIP_PROT = json.loads(r'''${skip_json}''')

    heavy_accs = { r["unimod"] for r in LABEL }
    if not heavy_accs:
        sys.exit("quantify_heavy_light: heavy_label is empty; nothing to pair.")

    base = pd.read_csv(BASE, sep='\\t')
    if base.empty:
        for p in ("peptide.abs_quant.long.tsv","protein.abs_quant.long.tsv",
                  "peptide.abs_quant.matrix.tsv","abs_quant_qc.tsv"):
            pd.DataFrame().to_csv(p, sep='\\t', index=False)
        sys.exit("no rows in base_long; empty outputs written.")

    heavy_row_count = int((base.get('channel', pd.Series(dtype=str)) == 'heavy').sum())
    if heavy_row_count == 0:
        sys.exit(
            "HARD ERROR: heavy_label is set but zero rows carry channel=='heavy' "
            "in base_quant.long.tsv. Diagnoses to check, in order:\\n"
            "  1. Heavy peptides are absent from the search library, and "
            "COMPLETE_CHANNELS did not run, or ran in the wrong direction.\\n"
            "  2. BLIB_TO_DLIB dropped a channel during conversion. Count "
            "entries either side: `SELECT COUNT(*) FROM RefSpectra` on the "
            ".complete.blib vs `SELECT COUNT(*) FROM entries` on the .dlib. "
            "A halving means BlibFile's inner join on "
            "SpectrumSourceFiles.id == RefSpectra.fileID discarded the rows "
            "COMPLETE_CHANNELS synthesised, which happens when their fileId "
            "is NULL (fixed in bin/complete_channels.py. A stale .complete.blib "
            "from before that fix will still show it; delete and regenerate). "
            "Note the drop hits the SYNTHESISED channel, so on a heavy-only "
            "vendor blib it removes the lights and this guard fires on the "
            "downstream ratio, not on the heavies themselves.\\n"
            "  3. A GPF-derived .elib was built on unspiked material. GPF "
            "injections must contain the spike-ins for absolute quant to work.\\n"
            "  4. The MERGE_QUANT channel detection regex missed the heavy "
            "form used by this engine (UniMod vs delta-mass).\\n"
            "  5. The acquisition schedule never targeted the heavy precursors.\\n"
            "Refusing to emit an empty absolute-quant table silently."
        )

    fallback = ['openswath','encyclopedia','diathem','tric']
    methods_present = [m for m in fallback
                       if f'abundance_{m}' in base.columns
                       and base[f'abundance_{m}'].notna().any()]
    if not methods_present:
        sys.exit("no abundance_<method> columns are populated in base_long.")

    _flank_re = re.compile(r'^[A-Z\\-]\\.|\\.[A-Z\\-]\$')
    def strip(s):
        s = re.sub(r'\\.\\([^)]+\\)|\\([^)]+\\)|\\[[^\\]]*\\]', '', str(s))
        return _flank_re.sub('', s).strip('.')
    if 'stripped_seq' not in base.columns:
        base['stripped_seq'] = base['peptide'].apply(strip)

    key = ['sample_id','stripped_seq','charge']

    def pair_for_method(method):
        col = f'abundance_{method}'
        if col not in base.columns:
            return None
        _extra = [c for c in ('id_qvalue','abundance_primary_source') if c in base.columns]
        sub = base[[*key,'channel','protein', col,'consistency', *_extra]].copy()

        if method == 'primary' and 'abundance_primary_source' in sub.columns:
            _src = (sub.dropna(subset=[col])
                       .groupby(key)['abundance_primary_source'].nunique())
            _bad = set(_src[_src > 1].index)
            if _bad:
                _mask = sub.set_index(key).index.isin(_bad)
                sub = sub[~_mask]
                print(f"pair_for_method(primary): dropped {len(_bad)} pair(s) whose two "
                      f"channels came from different methods", file=sys.stderr)
        if 'id_qvalue' not in sub.columns:
            sub['id_qvalue'] = pd.NA
        sub['id_qvalue'] = pd.to_numeric(sub['id_qvalue'], errors='coerce')

        if MIN_CONS > 0.0 and 'consistency' in sub.columns:
            sub.loc[sub['consistency'] < MIN_CONS, col] = pd.NA

        light = (sub[sub['channel']=='light']
                 .groupby(key).agg(light_amp=(col,'sum'),
                                    light_consistency=('consistency','mean'),
                                    light_qvalue=('id_qvalue','min'),
                                    protein=('protein','first'))
                 .reset_index())
        heavy = (sub[sub['channel']=='heavy']
                 .groupby(key).agg(heavy_amp=(col,'sum'),
                                    heavy_consistency=('consistency','mean'),
                                    heavy_qvalue=('id_qvalue','min'),
                                    heavy_protein=('protein','first'))
                 .reset_index())
        pr = light.merge(heavy, on=key, how='outer')
        pr['protein'] = pr['protein'].fillna(pr['heavy_protein'])
        pr = pr.drop(columns=['heavy_protein'])
        pr['method'] = method
        return pr

    paired_all = pd.concat([p for p in (pair_for_method(m) for m in methods_present + ['primary']) if p is not None],
                            ignore_index=True)

    std = pd.read_csv(STD, sep=None, engine='python')
    lc = {c.lower(): c for c in std.columns}

    def col_by(names):
        for n in names:
            if n in lc: return lc[n]
        return None

    seq_col  = col_by(('sequence','peptide','peptidesequence','compound','peptidoform'))
    chg_col  = col_by(('charge','z','precursorcharge'))
    amt_col  = col_by(('pmol','amount','heavy_pmol','known_pmol','amount per well [pmol]'))
    prot_col = col_by(('uniprot','protein','protein_id','accession','proteinaccession'))

    if amt_col is None:
        sys.exit(f"standard_amounts {STD}: need a pmol/amount column")

    std_ren = {amt_col: 'heavy_pmol_known'}

    if seq_col is not None:
        std['_seq_key'] = std[seq_col].apply(strip)
        std_ren[seq_col] = '_seq_key_orig'
        if chg_col is not None:
            std_ren[chg_col] = 'charge'
            cols = ['_seq_key','charge','heavy_pmol_known']
            paired_all = paired_all.merge(std.rename(columns=std_ren)[cols],
                left_on=['stripped_seq','charge'],
                right_on=['_seq_key','charge'], how='left')
        else:
            cols = ['_seq_key','heavy_pmol_known']
            paired_all = paired_all.merge(std.rename(columns=std_ren)[cols],
                left_on='stripped_seq', right_on='_seq_key', how='left')
    elif prot_col is not None:
        std_ren[prot_col] = '_prot_key'
        cols = ['_prot_key','heavy_pmol_known']
        prot_join = std.rename(columns=std_ren)[cols].dropna(subset=['_prot_key'])
        keys = list(prot_join['_prot_key'].astype(str))
        key2pmol = dict(zip(keys,
                             prot_join['heavy_pmol_known'].astype(float)))

        def find_pmol(prot_field):
            if pd.isna(prot_field): return None
            tokens = re.split(r'[,;|\\s]+', str(prot_field))
            for t in tokens:
                if t in key2pmol:
                    return key2pmol[t]
                head = t.split('_', 1)[0]
                if head in key2pmol:
                    return key2pmol[head]
            return None
        paired_all['heavy_pmol_known'] = paired_all['protein'].apply(find_pmol)
    else:
        sys.exit(f"standard_amounts {STD}: need either a sequence or a protein/uniprot column")

    _light = pd.to_numeric(paired_all['light_amp'], errors='coerce').replace(0.0, float('nan'))
    _heavy = pd.to_numeric(paired_all['heavy_amp'], errors='coerce').replace(0.0, float('nan'))
    paired_all['ratio_L_H']      = _light / _heavy
    paired_all['abs_pmol_light'] = paired_all['ratio_L_H'] * paired_all['heavy_pmol_known']

    def _absent(v):
        return pd.isna(v) or float(v) == 0.0
    def _skipped_standard(r):
        prot = str(r.get('protein') or '')
        return bool(prot) and any(p in prot for p in SKIP_PROT)
    def status(r):
        light_ok = not _absent(r['light_amp'])
        heavy_ok = not _absent(r['heavy_amp'])
        if light_ok and heavy_ok:
            return 'detected'
        if not light_ok and heavy_ok:
            return 'skipped_standard' if _skipped_standard(r) else 'light_missing'
        if light_ok and not heavy_ok:
            return 'heavy_missing'
        return 'both_missing'
    paired_all['status'] = paired_all.apply(status, axis=1)

    n_skipped = int((paired_all['status'] == 'skipped_standard').sum())
    if n_skipped:
        print(f"status skipped_standard: {n_skipped} rows matching "
              f"{SKIP_PROT} are heavy-only by design and are NOT counted "
              "as light_missing.", file=sys.stderr)

    n_abs = int(paired_all['abs_pmol_light'].notna().sum())
    if n_abs == 0:
        n_ratio = int(paired_all['ratio_L_H'].notna().sum())
        n_known = int(paired_all['heavy_pmol_known'].notna().sum())
        n_prot  = int(paired_all['protein'].astype(str).str.len().gt(0).sum()) \\
                  if 'protein' in paired_all.columns else 0
        sys.exit(
            "HARD ERROR: zero absolute amounts produced from "
            f"{len(paired_all)} paired rows.\\n"
            f"  light/heavy ratio computable : {n_ratio}\\n"
            f"  rows carrying a protein      : {n_prot}\\n"
            f"  rows matched to a known pmol : {n_known}\\n"
            "Read those three in order: the first zero is the broken link.\\n"
            "  * ratio 0        -> no peptide has BOTH channels quantified. "
            "Check abundance_<method> in base_quant.long.tsv; on the "
            "encyclopedia PRM path that comes from the .elib peptidequants "
            "table, which -libexport leaves near-empty when the search "
            "identified little.\\n"
            "  * protein 0      -> the peptide -> parent-sequence link is "
            "missing. With construct standards the pmol is "
            "priced per CONSTRUCT, so a peptide with no protein can never "
            "be costed. Check the protein column in base_quant.long.tsv and "
            "the peptidetoprotein table in analysis.elib.\\n"
            "  * known pmol 0   -> protein IDs are present but none matched "
            f"a key in {STD}. Compare the accession forms on both sides.\\n"
            "Refusing to emit an all-blank absolute-quant table silently."
        )
    print(f"absolute values produced: {n_abs} of {len(paired_all)} paired rows",
          file=sys.stderr)

    out_cols = ['sample_id','stripped_seq','charge','protein','method',
                'light_amp','heavy_amp','light_consistency','heavy_consistency',
                'ratio_L_H','heavy_pmol_known','abs_pmol_light','status']
    out_cols = [c for c in out_cols if c in paired_all.columns]
    paired_all[out_cols].to_csv("peptide.abs_quant.long.tsv", sep='\\t', index=False)

    pairs = paired_all[(paired_all['method'] == 'primary') &
                       (paired_all['status'] == 'detected')].copy()
    pairs['ratio_H_L'] = pairs['heavy_amp'] / pairs['light_amp'].replace(0, pd.NA)

    pairs['max_qvalue'] = pairs[['light_qvalue','heavy_qvalue']].max(axis=1)
    pairs['both_identified'] = pairs['max_qvalue'].le(FDR).fillna(False)
    pair_cols = ['sample_id','protein','stripped_seq','charge',
                 'light_amp','heavy_amp','ratio_L_H','ratio_H_L',
                 'heavy_pmol_known','abs_pmol_light',
                 'light_qvalue','heavy_qvalue','max_qvalue','both_identified',
                 'light_consistency','heavy_consistency']
    pair_cols = [c for c in pair_cols if c in pairs.columns]
    pairs = pairs[pair_cols].sort_values(['sample_id','protein','stripped_seq','charge'])
    pairs.to_csv("peptide.pairs.tsv", sep='\\t', index=False)
    conf = pairs[pairs['both_identified']]
    conf.to_csv("peptide.pairs.confident.tsv", sep='\\t', index=False)
    print(f"pairs table: {len(pairs)} light/heavy pairs with both channels quantified, "
          f"of {len(paired_all[paired_all['method']=='primary'])} primary-method rows",
          file=sys.stderr)
    print(f"  of those, BOTH CHANNELS IDENTIFIED at q<={FDR}: {len(conf)}", file=sys.stderr)
    for sid, g in conf.groupby('sample_id'):
        print(f"    {sid}: {len(g)} confident pairs, "
              f"{int(g['abs_pmol_light'].notna().sum())} with an absolute amount",
              file=sys.stderr)

    per_peptide_prot = (paired_all.dropna(subset=['abs_pmol_light'])
                                   .groupby(['sample_id','protein','method','stripped_seq'])
                                   .agg(abs_pmol_light=('abs_pmol_light','median'))
                                   .reset_index())
    prot = (per_peptide_prot.groupby(['sample_id','protein','method'])
                             .agg(n_peptides=('stripped_seq','nunique'),
                                  abs_pmol_light=('abs_pmol_light','median'))
                             .reset_index())
    prot.to_csv("protein.abs_quant.long.tsv", sep='\\t', index=False)

    primary_slice = paired_all[paired_all['method'] == 'primary']
    mat = primary_slice.pivot_table(index=['stripped_seq','charge','protein'],
                                     columns='sample_id',
                                     values='abs_pmol_light',
                                     aggfunc='first').reset_index()
    mat.to_csv("peptide.abs_quant.matrix.tsv", sep='\\t', index=False)

    qc = paired_all.dropna(subset=['heavy_pmol_known'])[
        ['sample_id','stripped_seq','charge','method',
         'heavy_pmol_known','abs_pmol_light','status']]
    qc.to_csv("abs_quant_qc.tsv", sep='\\t', index=False)

    print(f"methods evaluated: {methods_present + ['primary']}", file=sys.stderr)
    print(f"paired rows:       {len(paired_all)}", file=sys.stderr)
    for m in methods_present + ['primary']:
        sub = paired_all[paired_all['method']==m]
        n = int(sub['abs_pmol_light'].notna().sum())
        print(f"  {m:14}: {n} rows with abs value", file=sys.stderr)
    PY

    cat <<-EOF > versions.yml
    "${task.process}":
      python:  \$(python3 --version 2>&1 | awk '{print \$2}')
      primary: '${primary}'
      min_consistency: '${min_consistency}'
    EOF
    """

    stub:
    """
    touch peptide.pairs.tsv peptide.pairs.confident.tsv peptide.abs_quant.long.tsv protein.abs_quant.long.tsv \\
          peptide.abs_quant.matrix.tsv abs_quant_qc.tsv quantify_heavy_light.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
