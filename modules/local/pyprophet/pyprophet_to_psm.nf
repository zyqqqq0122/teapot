process PYPROPHET_TO_PSM {
    label 'python'
    tag   "${meta.id}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/pyprophet_prm/${meta.id}" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(reference_list), path(export_tsv)

    output:
    tuple val(meta), path("${meta.id}.psm.reference.txt"),     emit: psm
    tuple val(meta), path("${meta.id}.peptide.reference.txt"), emit: peptide
    tuple val(meta), path("${meta.id}.pyprophet_to_psm.log"),  emit: log
    path  "versions.yml", emit: versions

    script:
    """
    export HOME=\$PWD
    export MPLCONFIGDIR=\$PWD/.mplconfig
    export XDG_CACHE_HOME=\$PWD/.cache
    mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME"


    python3 - <<'PY' > ${meta.id}.pyprophet_to_psm.log 2>&1
    import pandas as pd, re, sys

    ML  = "${reference_list}"
    EXP = "${export_tsv}"
    SID = "${meta.id}"

    d = pd.read_csv(EXP, sep='\\t')
    print(f"pyprophet export rows: {len(d)}", file=sys.stderr)
    print(f"columns: {list(d.columns)}", file=sys.stderr)

    def pick(*names):
        for n in names:
            for c in d.columns:
                if c.lower() == n.lower():
                    return c
        return None

    pep_col  = pick('FullPeptideName', 'Sequence', 'peptide')
    prot_col = pick('ProteinName', 'Proteins', 'protein')
    q_col    = pick('m_score', 'q_value', 'qvalue')
    s_col    = pick('d_score', 'score')
    pep_prob = pick('pep', 'posterior_error_prob')
    dec_col  = pick('decoy')
    id_col   = pick('transition_group_id', 'SpecId', 'id')
    rt_col   = pick('RT', 'rt')
    if pep_col is None or q_col is None:
        sys.exit(f"PYPROPHET_TO_PSM: cannot find peptide/q-value columns in {EXP}")

    strip = lambda s: re.sub(r'\\(UniMod:\\d+\\)|\\[[^\\]]*\\]|\\.', '', str(s))

    ref = set()
    if ML and ML != 'NO_FILE':
        import csv as _csv
        with open(ML) as fh:
            sample = fh.read(4096); fh.seek(0)
            delim = '\\t' if sample.count('\\t') >= sample.count(',') else ','
            for row in _csv.reader(fh, delimiter=delim):
                if row and row[0] and not row[0].lower().startswith('compound'):
                    ref.add(strip(row[0]))
    print(f"reference peptides in mass list: {len(ref)}", file=sys.stderr)

    d['_seq'] = d[pep_col].apply(strip)
    if dec_col is not None:
        d['_decoy'] = pd.to_numeric(d[dec_col], errors='coerce').fillna(0).astype(int)
    else:
        d['_decoy'] = 0

    sub = d[d['_seq'].isin(ref)] if ref else d
    print(f"rows on the reference side: {len(sub)} "
          f"(targets {int((sub._decoy == 0).sum())}, decoys {int((sub._decoy == 1).sum())})",
          file=sys.stderr)

    out = pd.DataFrame({
        'SpecId': sub[id_col].astype(str) if id_col else [f"{SID}:{i}" for i in range(len(sub))],
        'ScanNr': pd.to_numeric(sub[rt_col], errors='coerce').fillna(0).astype(int) if rt_col else 0,
        'score':  pd.to_numeric(sub[s_col], errors='coerce') if s_col else 0.0,
        'q-value': pd.to_numeric(sub[q_col], errors='coerce'),
        'posterior_error_prob': (pd.to_numeric(sub[pep_prob], errors='coerce')
                                 if pep_prob else pd.NA),
        'Peptide': sub[pep_col].astype(str),
        'Proteins': sub[prot_col].astype(str) if prot_col else '',
    })
    out = out[sub['_decoy'].values == 0]
    out.to_csv(f"{SID}.psm.reference.txt", sep='\\t', index=False)

    pepout = (out.sort_values('q-value')
                 .drop_duplicates(subset=['Peptide'], keep='first'))
    pepout.to_csv(f"{SID}.peptide.reference.txt", sep='\\t', index=False)
    print(f"wrote {len(out)} PSM rows, {len(pepout)} peptide rows", file=sys.stderr)
    q = pd.to_numeric(out['q-value'], errors='coerce')
    print(f"passing q<=${params.fdr}: {int((q <= ${params.fdr}).sum())}", file=sys.stderr)
    PY

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
    EOF
    """

    stub:
    """
    touch ${meta.id}.psm.reference.txt ${meta.id}.peptide.reference.txt \\
          ${meta.id}.pyprophet_to_psm.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
