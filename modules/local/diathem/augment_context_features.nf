process AUGMENT_CONTEXT_FEATURES {
    label 'python'
    tag   "${meta.id}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/augment_context_features/${meta.id}" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(reference_tsv), path(background_tsv)
    path  diathem_tsv

    output:
    tuple val(meta), path("${meta.id}_reference.augmented.tsv"),
                     path("${meta.id}_background.augmented.tsv"), emit: split
    path  "${meta.id}.augment.log",                                emit: log
    path  "versions.yml", emit: versions

    script:
    """

    python3 - <<'PY' > ${meta.id}.augment.log 2>&1
    import pandas as pd, os, re, sys

    REF  = "${reference_tsv}"
    BG   = "${background_tsv}"
    DIA  = "${diathem_tsv}"
    SID  = "${meta.id}"
    OUT_REF = "${meta.id}_reference.augmented.tsv"
    OUT_BG  = "${meta.id}_background.augmented.tsv"

    def strip_mods(s):
        return re.sub(r'\\.\\([^)]+\\)|\\([^)]+\\)|\\[[^\\]]*\\]', '', str(s)).strip('.')

    if os.path.basename(DIA) == 'NO_FILE' or not os.path.getsize(DIA):
        dia = pd.DataFrame()
    else:
        dia = pd.read_csv(DIA, sep='\\t')
        dia = dia.rename(columns={'run_id': 'sample_id'})
        dia = dia[dia['sample_id'] == SID].copy()
        dia['stripped_seq'] = dia['peptidoform'].apply(strip_mods)
        dia = dia[['stripped_seq','charge','consistency','n_effective_transitions']]
        dia = (dia.groupby(['stripped_seq','charge'], as_index=False)
                  .agg(consistency=('consistency','median'),
                       n_effective_transitions=('n_effective_transitions','median')))

    def augment(path, out_path):
        t = pd.read_csv(path, sep='\\t')
        if t.empty:
            t.to_csv(out_path, sep='\\t', index=False); return 0, 0
        pep_col_idx = len(t.columns) - 2

        t['_seq_key'] = t['Peptide'].apply(strip_mods)
        def pep_to_charge(p):
            m = re.search(r'(?:\\+|/)(\\d+)\$', str(p))
            return int(m.group(1)) if m else 2
        t['_charge'] = t['Peptide'].apply(pep_to_charge)

        if dia.empty:
            t['xrunConsist']    = 0.0
            t['xrunEffTrans']   = 0.0
            t['xrunNA']         = 1
            joined = 0
        else:
            m = t.merge(dia, left_on=['_seq_key','_charge'],
                             right_on=['stripped_seq','charge'],
                             how='left', validate='many_to_one')
            joined = int(m['consistency'].notna().sum())
            med_c  = m['consistency'].median() if joined else 0.0
            med_e  = m['n_effective_transitions'].median() if joined else 0.0
            t['xrunConsist']  = m['consistency'].fillna(med_c).values
            t['xrunEffTrans'] = m['n_effective_transitions'].fillna(med_e).values
            t['xrunNA']       = m['consistency'].isna().astype(int).values

        t = t.drop(columns=['_seq_key','_charge'])

        new_cols = ['xrunConsist','xrunEffTrans','xrunNA']
        base_cols = [c for c in t.columns if c not in new_cols]
        pep_col = base_cols[-2]
        prot_col = base_cols[-1]
        reorder = base_cols[:-2] + new_cols + [pep_col, prot_col]
        t[reorder].to_csv(out_path, sep='\\t', index=False)
        return len(t), joined

    n_ref, j_ref = augment(REF, OUT_REF)
    n_bg,  j_bg  = augment(BG,  OUT_BG)
    print(f"reference: {n_ref} rows, {j_ref} joined to diathem", file=sys.stderr)
    print(f"background: {n_bg} rows, {j_bg} joined to diathem", file=sys.stderr)
    PY

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
    EOF
    """

    stub:
    """
    touch ${meta.id}_reference.augmented.tsv \\
          ${meta.id}_background.augmented.tsv \\
          ${meta.id}.augment.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
