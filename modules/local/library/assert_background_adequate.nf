process ASSERT_BACKGROUND_ADEQUATE {
    label 'python'
    tag   "${union_dlib.baseName}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/merge_libraries" },
               mode: params.publish_mode

    input:
    path union_dlib
    path reference_list
    val  min_targets

    output:
    path "background_adequacy.tsv", emit: report
    path "background_adequacy.log", emit: log
    path "versions.yml",            emit: versions

    script:
    """
    python3 - <<'PY' > background_adequacy.log 2>&1
    import sqlite3, csv, re, sys

    UNION_DLIB      = "${union_dlib}"
    REFERENCE_LIST  = "${reference_list}"
    MIN_TARGETS     = int("${min_targets}")

    def strip_mods(s):
        return re.sub(r'\\[[^\\]]*\\]', '', str(s)).strip()

    with open(REFERENCE_LIST) as f:
        head = f.readline()
    delim = '\\t' if head.count('\\t') > head.count(',') else ','
    with open(REFERENCE_LIST) as f:
        r = csv.reader(f, delimiter=delim)
        header = next(r)
        seq_idx = None
        for i, c in enumerate(header):
            if c.strip().lower() in ('compound','peptidesequence','sequence','peptide'):
                seq_idx = i; break
        if seq_idx is None:
            sys.exit(f"reference_list {REFERENCE_LIST}: no Compound/PeptideSequence column")
        ref_seqs = set()
        for row in r:
            if len(row) <= seq_idx: continue
            s = strip_mods(row[seq_idx])
            if s:
                ref_seqs.add(s)
    n_ref_peptides = len(ref_seqs)

    con = sqlite3.connect(UNION_DLIB)
    union_seqs = set()
    for (seq,) in con.execute("SELECT DISTINCT PeptideSeq FROM entries"):
        union_seqs.add(seq)
    con.close()
    n_union_peptides = len(union_seqs)

    bg_seqs = union_seqs - ref_seqs
    n_bg_peptides = len(bg_seqs)

    with open("background_adequacy.tsv", "w", newline='') as f:
        w = csv.writer(f, delimiter='\\t')
        w.writerow(["metric", "value"])
        w.writerow(["reference_peptides", n_ref_peptides])
        w.writerow(["union_peptides",     n_union_peptides])
        w.writerow(["background_peptides", n_bg_peptides])
        w.writerow(["min_background_required", MIN_TARGETS])

    print(f"reference peptides:         {n_ref_peptides}", file=sys.stderr)
    print(f"union peptides:             {n_union_peptides}", file=sys.stderr)
    print(f"background peptides:        {n_bg_peptides}", file=sys.stderr)
    print(f"min background required:    {MIN_TARGETS}", file=sys.stderr)

    if n_bg_peptides < MIN_TARGETS:
        sys.exit(
            f"HARD ERROR: background adequacy check failed. "
            f"union has {n_union_peptides} peptides, reference has "
            f"{n_ref_peptides}, so background = {n_bg_peptides}, which is"
            f"below the required {MIN_TARGETS}. Context's mProphet LDA "
            f"trains on background target/decoy separation and needs "
            f"a meaningful number of background targets. Likely causes:\\n"
            f"  * background_source is 'none' or an unset library;\\n"
            f"  * background_source='koina' but background_fasta is the "
            f"vendor construct FASTA (predicting from vendor gives all "
            f"vendor targets, which collide with the target library);\\n"
            f"  * the background library happens to be a subset of the "
            f"target library.\\n"
            f"Tune background_min_targets down only if you know the "
            f"downstream Context model will still be non-degenerate."
        )
    PY

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
    EOF
    """

    stub:
    """
    touch background_adequacy.tsv background_adequacy.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
