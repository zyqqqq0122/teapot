process SPLIT_OPENSWATH_FEATURES {
    label 'python'
    tag   "${meta.id}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/split_openswath_features/${meta.id}" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(reference_list), path(osw)

    output:
    tuple val(meta), path("${meta.id}_reference.tsv"),
                     path("${meta.id}_background.tsv"),     emit: split
    tuple val(meta), path("${meta.id}.split.log"),          emit: log
    path  "versions.yml", emit: versions

    script:
    """

    python3 - <<'PY' > ${meta.id}.split.log 2>&1
    import sqlite3, re, sys, csv
    import pandas as pd

    OSW  = "${osw}"
    ML   = "${reference_list}"
    REF_OUT = "${meta.id}_reference.tsv"
    BG_OUT  = "${meta.id}_background.tsv"

    def strip_mods(s):
        return re.sub(r'\\[[^\\]]*\\]', '', str(s)).strip()

    with open(ML) as f:
        head = f.readline()
    delim = '\\t' if head.count('\\t') > head.count(',') else ','
    ml_df = pd.read_csv(ML, sep=delim)
    seq_col = None
    for c in ml_df.columns:
        if c.strip().lower() in ('compound', 'peptidesequence', 'sequence', 'peptide'):
            seq_col = c; break
    if seq_col is None:
        sys.exit(f"mass list {ML}: no Compound/PeptideSequence column found")
    ref_seqs = { strip_mods(x) for x in ml_df[seq_col].dropna() }
    print(f"reference peptides in mass list: {len(ref_seqs)}", file=sys.stderr)

    con = sqlite3.connect(OSW)

    def var_cols(table):
        cur = con.execute(f"PRAGMA table_info({table})")
        return [r[1] for r in cur if r[1].lower().startswith('var_')]

    ms2_vars = var_cols('FEATURE_MS2')
    try:
        ms1_vars = var_cols('FEATURE_MS1')
    except sqlite3.OperationalError:
        ms1_vars = []

    ms2_sel = ','.join(f'FMS2.{c} AS {c}' for c in ms2_vars)
    ms1_sel = ','.join(f'FMS1.{c} AS ms1_{c}' for c in ms1_vars)

    joins = "JOIN FEATURE_MS2 FMS2 ON FMS2.FEATURE_ID = F.ID"
    if ms1_vars:
        joins += " LEFT JOIN FEATURE_MS1 FMS1 ON FMS1.FEATURE_ID = F.ID"

    select_parts = [
        "F.ID AS feature_id",
        "F.RUN_ID AS run_id",
        "F.EXP_RT AS exp_rt",
        "P.ID AS peptide_id",
        "P.UNMODIFIED_SEQUENCE AS peptide_seq",
        "P.MODIFIED_SEQUENCE AS peptide_modseq",
        "P.DECOY AS decoy",
        "PRE.CHARGE AS charge",
        "PRE.TRAML_ID AS traml_id",
    ]
    if ms1_sel:
        select_parts.append(ms1_sel)
    select_parts.append(ms2_sel)

    sql = f'''
        SELECT {",".join(select_parts)}
        FROM FEATURE F
        JOIN PRECURSOR_PEPTIDE_MAPPING PPM ON PPM.PRECURSOR_ID = F.PRECURSOR_ID
        JOIN PEPTIDE P                     ON P.ID              = PPM.PEPTIDE_ID
        JOIN PRECURSOR PRE                 ON PRE.ID            = F.PRECURSOR_ID
        {joins}
    '''
    df = pd.read_sql_query(sql, con)
    print(f"features loaded: {len(df)}", file=sys.stderr)

    try:
        prot = pd.read_sql_query(
            "SELECT PPM.PEPTIDE_ID AS peptide_id, PROT.PROTEIN_ACCESSION AS accession "
            "FROM PEPTIDE_PROTEIN_MAPPING PPM "
            "JOIN PROTEIN PROT ON PROT.ID = PPM.PROTEIN_ID", con)
        prot_map = (prot.groupby('peptide_id')['accession']
                    .apply(lambda xs: ','.join(sorted(set(xs))))
                    .to_dict())
    except Exception:
        prot_map = {}
    con.close()

    con2 = sqlite3.connect(OSW)

    df['proteins'] = df['peptide_id'].map(prot_map).fillna(df['peptide_seq'])
    df['label']   = df['decoy'].apply(lambda d: -1 if int(d) == 1 else 1)
    df['SpecId']  = ('f' + df['feature_id'].astype(str) + '_' +
                     df['peptide_modseq'].astype(str) + '_' +
                     df['charge'].astype(str))
    df['ScanNr']  = df['exp_rt'].round().astype('Int64')
    df['Peptide'] = df['peptide_modseq']
    df['Proteins'] = df['proteins']

    target_map = pd.read_sql_query(
        "SELECT PRE.TRAML_ID AS traml_id, P.UNMODIFIED_SEQUENCE AS seq "
        "FROM PRECURSOR PRE "
        "JOIN PRECURSOR_PEPTIDE_MAPPING M ON M.PRECURSOR_ID = PRE.ID "
        "JOIN PEPTIDE P                   ON P.ID           = M.PEPTIDE_ID "
        "WHERE PRE.DECOY = 0", con2)
    traml_is_ref = dict(zip(target_map['traml_id'].astype(str),
                            target_map['seq'].isin(ref_seqs)))
    con2.close()

    is_decoy   = df['decoy'].astype(int) == 1
    paired_id  = df['traml_id'].astype(str).str.replace('^DECOY_', '', regex=True)
    decoy_ref  = paired_id.map(traml_is_ref).fillna(False).astype(bool)
    target_ref = df['peptide_seq'].isin(ref_seqs)
    df['in_ref'] = target_ref.where(~is_decoy, decoy_ref)

    n_dec = int(is_decoy.sum())
    n_dec_paired = int((is_decoy & df['traml_id'].astype(str).str.startswith('DECOY_')).sum())
    if n_dec and n_dec_paired == 0:
        print("WARNING: no decoy precursor uses the 'DECOY_<target TRAML_ID>' "
              "naming this pairing relies on. The library was probably not built "
              "by OpenSwathDecoyGenerator (e.g. a user-supplied --openswath_pqp). "
              "Reference decoys cannot be identified, so Context's reference "
              "q-values will be uninformative (all 1).", file=sys.stderr)

    feat_cols = [c for c in df.columns
                 if c.lower().startswith(('var_', 'ms1_var_'))]

    df = df.rename(columns={'label': 'Label'})
    ordered = ['SpecId', 'Label', 'ScanNr'] + feat_cols + ['Peptide', 'Proteins']

    ref_df = df[df['in_ref']][ordered]
    bg_df  = df[~df['in_ref']][ordered]

    def counts(x, tag):
        t = int((x['Label'] ==  1).sum())
        d = int((x['Label'] == -1).sum())
        print(f"{tag}: rows={len(x)} targets={t} decoys={d}", file=sys.stderr)
        return t, d

    ref_t, ref_d = counts(ref_df, "reference")
    bg_t,  bg_d  = counts(bg_df,  "background")

    if ref_d == 0:
        print("WARNING: the reference set contains ZERO decoys. Context cannot "
              "estimate a reference FDR without them and will return q-value = 1 "
              "for every reference peptide, regardless of how well they score. "
              "Any downstream abundance is therefore unvalidated.", file=sys.stderr)

    ref_df.to_csv(REF_OUT, sep='\\t', index=False)
    bg_df.to_csv (BG_OUT,  sep='\\t', index=False)

    if bg_t == 0:
        sys.exit(
            f"HARD ERROR: SPLIT_OPENSWATH_FEATURES produced zero background targets "
            f"for {REF_OUT.split('_reference')[0]}. The reference list appears "
            f"to cover the entire search library, so Context's LDA has nothing "
            f"to train on. Fix by unioning a background library into the search "
            f"transitions. For the OSW blib route set background_koina=true "
            f"(a Koina-predicted proteome library will be unioned with the "
            f"vendor blib); reference_list is derived from the vendor blib only."
        )
    PY

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
    EOF
    """

    stub:
    """
    touch ${meta.id}_reference.tsv ${meta.id}_background.tsv ${meta.id}.split.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
