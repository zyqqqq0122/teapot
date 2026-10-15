process QUANTIFY_HEAVY_LIGHT {
    label 'python'
    tag   "abs_quant"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/quant/${route}" },
               mode: params.publish_mode,
               saveAs: { fn -> fn == 'versions.yml' ? 'quantify_heavy_light.versions.yml' : fn }

    input:
    path  base_long
    path  standard_amounts
    path  fasta, stageAs: 'fasta/*'
    val   heavy_label
    val   primary
    val   min_consistency
    val   route

    output:
    path "peptide.abs_quant.tsv",           emit: peptide
    path "peptide.abs_quant.confident.tsv", emit: peptide_confident
    path "peptide.abs_quant.matrix.tsv",    emit: peptide_matrix
    path "protein.abs_quant.tsv",           emit: protein
    path "quantify_heavy_light.log",        emit: log
    path "versions.yml",                    emit: versions

    script:
    def label_json = groovy.json.JsonOutput.toJson(heavy_label ?: [])
    def skip_json  = groovy.json.JsonOutput.toJson(params.complete_channels_skip_proteins ?: [])
    """

    python3 - <<'PY' > quantify_heavy_light.log 2>&1
    import pandas as pd, re, sys, json

    BASE      = "${base_long}"
    STD       = "${standard_amounts}"
    FASTA     = "${fasta}"
    LABEL     = json.loads(r'''${label_json}''')
    PRIMARY   = "${primary}"
    MIN_CONS  = float("${min_consistency}")
    FDR       = float("${params.fdr}")
    SKIP_PROT = json.loads(r'''${skip_json}''')
    TOL       = float("${params.pair_rt_tolerance_seconds ?: 0}")

    heavy_accs = { r["unimod"] for r in LABEL }
    if not heavy_accs:
        sys.exit("quantify_heavy_light: heavy_label is empty; nothing to pair.")

    base = pd.read_csv(BASE, sep='\\t')
    if base.empty:
        sys.exit("quantify_heavy_light: %s has no rows." % BASE)

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

    methods = [m for m in ('encyclopedia', 'openswath', 'diathem', 'tric')
               if f'abundance_{m}' in base.columns and base[f'abundance_{m}'].notna().any()]
    if not methods:
        sys.exit("no abundance_<method> columns are populated in base_long.")
    if PRIMARY not in methods:
        sys.exit(f"primary_abundance = '{PRIMARY}' has no abundances on this route; "
                 f"tools with abundances: {methods}.")

    _flank_re = re.compile(r'^[A-Z\\-]\\.|\\.[A-Z\\-]\$')
    def strip(s):
        s = _flank_re.sub('', str(s).strip())
        s = re.sub(r'\\.\\([^)]+\\)|\\([^)]+\\)|\\[[^\\]]*\\]', '', s)
        return s.strip('.')
    base['peptide'] = base['peptidoform'].fillna(base['peptide']).apply(strip)

    def in_fasta(path, peptides):
        peptides = {p for p in peptides if p}
        short = {p for p in peptides if len(p) < 5}
        by_start = {}
        for p in peptides - short:
            by_start.setdefault(p[:5], []).append(p)
        found = {}

        def scan(name, s):
            for i in range(len(s) - 4):
                for p in by_start.get(s[i:i + 5], ()):
                    if s.startswith(p, i):
                        found.setdefault(p, set()).add(name)
            for p in short:
                if p in s:
                    found.setdefault(p, set()).add(name)
        name, seq = None, []
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if name is not None:
                        scan(name, ''.join(seq).upper())
                    name, seq = line[1:].split(' ', 1)[0], []
                else:
                    seq.append(line)
        if name is not None:
            scan(name, ''.join(seq).upper())
        return {p: ';'.join(sorted(n)) for p, n in found.items()}

    base['protein'] = (base['protein'] if 'protein' in base.columns else pd.NA)
    base['protein'] = base['protein'].astype(object)
    blank = base['protein'].isna() | base['protein'].astype(str).str.strip().eq('')
    if blank.any():
        found = in_fasta(FASTA, set(base.loc[blank, 'peptide']))
        base.loc[blank, 'protein'] = base.loc[blank, 'peptide'].map(found)
        print(f"protein: from the FASTA for {len(found)} of {base.loc[blank, 'peptide'].nunique()} "
              f"peptide(s) without one", file=sys.stderr)
    for c in ('consistency', 'rt_apex_seconds'):
        if c not in base.columns:
            base[c] = pd.NA
    base['id_qvalue'] = pd.to_numeric(base.get('id_qvalue'), errors='coerce')
    if 'peak_group_rank' in base.columns:
        rank = pd.to_numeric(base['peak_group_rank'], errors='coerce')
        base = base[rank.isna() | (rank <= 1)]

    pc = ['peptide', 'charge']
    has_std = pd.MultiIndex.from_frame(base[pc]).isin(
        pd.MultiIndex.from_frame(base.loc[base['channel'] == 'heavy', pc]))
    if not has_std.all():
        print(f"left out {base.loc[~has_std, pc].drop_duplicates().shape[0]} precursor(s) "
              f"with no heavy row in any run: no standard", file=sys.stderr)
    base = base[has_std]

    key = ['sample_id', 'peptide', 'charge']

    def rt_col(method):
        c = f'rt_{method}_seconds'
        return c if c in base.columns and base[c].notna().any() else 'rt_apex_seconds'

    seen = [c for c in ('id_qvalue', 'id_pep', 'rt_apex_seconds', 'rt_encyclopedia_seconds',
                        'abundance_encyclopedia', 'abundance_openswath', 'abundance_tric')
            if c in base.columns]

    def pair_for_method(method):
        col = f'abundance_{method}'
        rows = base if method == 'diathem' else base[base[seen].notna().any(axis=1)]
        sub = rows[[*key, 'channel', 'protein', col, 'consistency', 'id_qvalue']].copy()
        sub['rt'] = pd.to_numeric(rows[rt_col(method)], errors='coerce')
        sub[col] = pd.to_numeric(sub[col], errors='coerce')
        if f'{col}_shared_ions' in rows.columns:
            shared = pd.to_numeric(rows[f'{col}_shared_ions'], errors='coerce')
            sub[col] = shared.where(shared.notna(), sub[col])
            print(f"{method}: where both forms are quantified ({int(shared.notna().sum())} rows), "
                  f"each is summed over the quant ions light and heavy share", file=sys.stderr)
        if MIN_CONS > 0.0:
            sub.loc[pd.to_numeric(sub['consistency'], errors='coerce') < MIN_CONS, col] = pd.NA
        agg = lambda pre: dict(**{f'{pre}_abundance': (col, lambda s: s.sum(min_count=1)),
                                  f'{pre}_q_value': ('id_qvalue', 'min'),
                                  f'{pre}_consistency': ('consistency', 'mean'),
                                  f'{pre}_protein': ('protein', 'first'),
                                  f'{pre}_rt_seconds': ('rt', 'first')})
        light = sub[sub['channel'] == 'light'].groupby(key).agg(**agg('light')).reset_index()
        heavy = sub[sub['channel'] == 'heavy'].groupby(key).agg(**agg('heavy')).reset_index()
        pr = light.merge(heavy, on=key, how='outer')
        pr['protein'] = pr['light_protein'].fillna(pr['heavy_protein'])
        pr['method'] = method
        return pr.drop(columns=['light_protein', 'heavy_protein'])

    pairs = pd.concat([pair_for_method(m) for m in methods], ignore_index=True)

    if SKIP_PROT:
        skipped = pairs['protein'].astype(str).apply(lambda p: any(s in p for s in SKIP_PROT))
        n_skip = int(skipped.sum())
        pairs = pairs[~skipped].copy()
        if n_skip:
            print(f"left out {n_skip} pair rows of {SKIP_PROT}, heavy-only by design",
                  file=sys.stderr)

    gap = (pairs['light_rt_seconds'] - pairs['heavy_rt_seconds']).abs()
    pairs['co_eluting'] = pd.array([pd.NA if (TOL <= 0 or pd.isna(g)) else bool(g <= TOL)
                                    for g in gap], dtype='boolean')

    std = pd.read_csv(STD, sep=None, engine='python')
    lc = {c.lower(): c for c in std.columns}
    def col_by(names):
        return next((lc[n] for n in names if n in lc), None)
    seq_col  = col_by(('sequence', 'peptide', 'peptidesequence', 'compound', 'peptidoform'))
    chg_col  = col_by(('charge', 'z', 'precursorcharge'))
    amt_col  = col_by(('pmol', 'amount', 'heavy_pmol', 'known_pmol', 'amount per well [pmol]'))
    prot_col = col_by(('uniprot', 'protein', 'protein_id', 'accession', 'proteinaccession'))
    if amt_col is None:
        sys.exit(f"standard_amounts {STD}: need a pmol/amount column")

    if seq_col is not None:
        known = pd.DataFrame({'peptide': std[seq_col].apply(strip),
                              'heavy_pmol_known': pd.to_numeric(std[amt_col], errors='coerce')})
        on = ['peptide']
        if chg_col is not None:
            known['charge'] = pd.to_numeric(std[chg_col], errors='coerce')
            on.append('charge')
        pairs = pairs.merge(known.drop_duplicates(subset=on), on=on, how='left')
    elif prot_col is not None:
        key2pmol = dict(zip(std[prot_col].astype(str),
                            pd.to_numeric(std[amt_col], errors='coerce')))
        def find_pmol(prot_field):
            if pd.isna(prot_field): return None
            for t in re.split(r'[,;|\\s]+', str(prot_field)):
                if t in key2pmol: return key2pmol[t]
                head = t.split('_', 1)[0]
                if head in key2pmol: return key2pmol[head]
            return None
        pairs['heavy_pmol_known'] = pairs['protein'].apply(find_pmol)
    else:
        sys.exit(f"standard_amounts {STD}: need either a sequence or a protein/uniprot column")

    light = pairs['light_abundance'].replace(0.0, float('nan'))
    heavy = pairs['heavy_abundance'].replace(0.0, float('nan'))
    pairs['ratio_L_H'] = light / heavy
    pairs['ratio_H_L'] = heavy / light
    pairs['abs_pmol_light'] = pairs['ratio_L_H'] * pd.to_numeric(pairs['heavy_pmol_known'], errors='coerce')
    pairs['light_detected'] = pairs['light_q_value'] <= FDR
    pairs['heavy_detected'] = pairs['heavy_q_value'] <= FDR

    def status(r):
        l, h = pd.notna(light[r.name]), pd.notna(heavy[r.name])
        if l and h: return 'both_quantified'
        if h: return 'light_missing'
        return 'heavy_missing' if l else 'both_missing'
    pairs['status'] = pairs.apply(status, axis=1)
    both_detected = pairs['light_detected'] & pairs['heavy_detected']
    apart = pairs['co_eluting'].eq(False).fillna(False).astype(bool)
    pairs['confident'] = (pairs['method'] == PRIMARY) & (pairs['status'] == 'both_quantified') & \\
                         both_detected & ~apart

    n_abs = int(pairs['abs_pmol_light'].notna().sum())
    if n_abs == 0:
        n_ratio = int(pairs['ratio_L_H'].notna().sum())
        n_known = int(pairs['heavy_pmol_known'].notna().sum())
        n_prot  = int(pairs['protein'].astype(str).str.len().gt(0).sum())
        sys.exit(
            "HARD ERROR: zero absolute amounts produced from "
            f"{len(pairs)} light/heavy pairs.\\n"
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

    cols = ['sample_id', 'method', 'protein', 'peptide', 'charge', 'status', 'confident',
            'light_abundance', 'heavy_abundance', 'ratio_L_H', 'ratio_H_L',
            'heavy_pmol_known', 'abs_pmol_light',
            'light_q_value', 'heavy_q_value', 'light_detected', 'heavy_detected',
            'light_rt_seconds', 'heavy_rt_seconds', 'co_eluting',
            'light_consistency', 'heavy_consistency']
    if 'diathem' not in methods:
        cols = [c for c in cols if not c.endswith('_consistency')]
    pairs = pairs[cols].sort_values(['method', 'sample_id', 'protein', 'peptide', 'charge'])
    pairs.to_csv("peptide.abs_quant.tsv", sep='\\t', index=False)

    conf = pairs[pairs['confident']].drop(columns=['method', 'status', 'confident', 'light_detected',
                                                   'heavy_detected', 'co_eluting'])
    conf.to_csv("peptide.abs_quant.confident.tsv", sep='\\t', index=False)

    amounts = conf.dropna(subset=['abs_pmol_light'])
    amounts.pivot_table(index=['protein', 'peptide', 'charge'], columns='sample_id',
                        values='abs_pmol_light', aggfunc='first') \\
        .reset_index().to_csv("peptide.abs_quant.matrix.tsv", sep='\\t', index=False)

    per_peptide = (amounts.groupby(['sample_id', 'protein', 'peptide'])
                       .agg(abs_pmol_light=('abs_pmol_light', 'median')).reset_index())
    (per_peptide.groupby(['sample_id', 'protein'])
                .agg(n_peptides=('peptide', 'nunique'),
                     abs_pmol_light=('abs_pmol_light', 'median'))
                .reset_index()
                .to_csv("protein.abs_quant.tsv", sep='\\t', index=False))

    print(f"light/heavy pairs per tool (q <= {FDR} on both channels = detected):", file=sys.stderr)
    for m, g in pairs.groupby('method'):
        q = g[g['status'] == 'both_quantified']
        d = q[q['light_detected'] & q['heavy_detected']]
        mark = '  <- primary_abundance: confident table, matrix, protein amounts' if m == PRIMARY else ''
        print(f"  {m:14} {q['peptide'].nunique():4} peptides with both channels quantified, "
              f"{d['peptide'].nunique():4} with both detected{mark}", file=sys.stderr)
        if TOL > 0:
            print(f"  {'':14} detected pair rows: {len(d)}; light and heavy peaks within "
                  f"{TOL:g} s: {int(d['co_eluting'].eq(True).sum())}, further apart: "
                  f"{int(d['co_eluting'].eq(False).sum())}, no retention time: "
                  f"{int(d['co_eluting'].isna().sum())}", file=sys.stderr)
    if TOL <= 0:
        print("co-elution of light and heavy not checked (pair_rt_tolerance_seconds 0)",
              file=sys.stderr)
    print(f"confident pairs: {len(conf)}, {len(amounts)} of them with a known amount "
          f"(matrix and protein amounts)", file=sys.stderr)
    if conf.empty:
        print(f"WARNING: no pair of '{PRIMARY}' has both channels detected at one time, so the "
              f"confident table, matrix and protein amounts are empty; "
              f"peptide.abs_quant.tsv has every tool's pairs.", file=sys.stderr)
    PY

    cat <<-EOF > versions.yml
    "${task.process}":
      python:  \$(python3 --version 2>&1 | awk '{print \$2}')
      primary: '${primary}'
      min_consistency: '${min_consistency}'
      pair_rt_tolerance_seconds: '${params.pair_rt_tolerance_seconds}'
    EOF
    """

    stub:
    """
    touch peptide.abs_quant.tsv peptide.abs_quant.confident.tsv peptide.abs_quant.matrix.tsv \\
          protein.abs_quant.tsv quantify_heavy_light.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
