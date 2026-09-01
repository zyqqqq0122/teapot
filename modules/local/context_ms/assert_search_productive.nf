process ASSERT_SEARCH_PRODUCTIVE {
    label 'python'
    tag   "${meta.id}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/qc/${meta.id}" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(reference_report)
    val   engine

    output:
    tuple val(meta), path("${meta.id}.${engine}.search_qc.tsv"), emit: qc
    path  "versions.yml", emit: versions

    script:
    """
    python3 - <<'PY' > ${meta.id}.${engine}.search_qc.tsv
    import csv, sys

    REPORT = "${reference_report}"
    FDR    = float("${params.fdr}")
    MIN_ID = int("${params.search_min_reference_ids}")
    ENGINE = "${engine}"
    SID    = "${meta.id}"

    with open(REPORT) as f:
        rows = [r for r in csv.DictReader(f, delimiter='\\t')]

    qcol = next((c for c in (rows[0] if rows else {})
                 if c.strip().lower() in ('q-value', 'q_value', 'qvalue')), None)
    n = len(rows)
    qs = []
    if qcol:
        for r in rows:
            v = (r.get(qcol) or '').strip()
            if v and v.upper() not in ('NA', 'NAN'):
                try:
                    qs.append(float(v))
                except ValueError:
                    pass

    passing = sum(1 for q in qs if q <= FDR)
    all_one = bool(qs) and all(q >= 1.0 for q in qs)

    w = csv.writer(sys.stdout, delimiter='\\t')
    w.writerow(['sample_id', 'engine', 'reference_rows', 'q_values_present',
                'passing_at_fdr', 'fdr', 'all_q_equal_1'])
    w.writerow([SID, ENGINE, n, len(qs), passing, FDR, all_one])
    sys.stdout.flush()

    if not rows:
        sys.exit("ASSERT_SEARCH_PRODUCTIVE [%s/%s]: the reference report is empty." % (SID, ENGINE))
    if qcol is None:
        sys.exit("ASSERT_SEARCH_PRODUCTIVE [%s/%s]: no q-value column in %s; cannot "
                 "verify the search produced identifications." % (SID, ENGINE, REPORT))
    if all_one:
        sys.exit(
            "ASSERT_SEARCH_PRODUCTIVE [%s/%s]: every reference q-value is exactly 1 "
            "over %d rows, so no FDR was estimated and nothing downstream is "
            "validated.\\n"
            "  Almost always this means the reference set has ZERO DECOYS -- with no "
            "decoys the q-value is 1 by construction whatever the scores do.\\n"
            "  encyclopedia: the reference list needs decoys "
            "(add_decoys_to_reference_list=true, or supply a mass list with an isDecoy column).\\n"
            "  openswath: check the SPLIT_OPENSWATH_FEATURES log for the reference decoy "
            "count; decoys are paired to targets via PRECURSOR.TRAML_ID = "
            "'DECOY_'+<target id>, which needs a library built by OpenSwathDecoyGenerator."
            % (SID, ENGINE, n))
    if passing == 0:
        sys.exit(
            "ASSERT_SEARCH_PRODUCTIVE [%s/%s]: ZERO reference peptides at q<=%g, out "
            "of %d scored. The search ran but identified nothing, so any abundance "
            "downstream would be unvalidated noise.\\n"
            "  First thing to check is INSTRUMENT-APPROPRIATE SEARCH PARAMETERS. "
            "EncyclopeDIA defaults to -ftol 10 -ftolunits ppm and OpenSWATH to "
            "-mz_extraction_window 30 ppm; both are Orbitrap settings. On ion-trap "
            "MS2 (mzML filter strings like 'ITMS + c NSI ... ms2 <mz>@hcd30') the "
            "fragment error is ~0.1-0.3 Da, the match window is too narrow, "
            "and the search returns nothing.\\n"
            "  encyclopedia: set context_search_args, e.g. "
            "'-ftol 0.4 -ftolunits AMU -lftol 0.4 -lftolunits AMU -frag HCD'.\\n"
            "  openswath: set -mz_extraction_window / -irt_mz_extraction_window in Th."
            % (SID, ENGINE, FDR, n))
    if passing < MIN_ID:
        print("WARNING: ASSERT_SEARCH_PRODUCTIVE [%s/%s]: only %d reference peptides "
              "at q<=%g (threshold %d). Legitimate for a very small assay, otherwise "
              "check the search parameters." % (SID, ENGINE, passing, FDR, MIN_ID),
              file=sys.stderr)
    PY

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
      fdr: ${params.fdr}
    EOF
    """

    stub:
    """
    echo -e "sample_id\\tengine\\treference_rows\\tq_values_present\\tpassing_at_fdr\\tfdr\\tall_q_equal_1" > ${meta.id}.${engine}.search_qc.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
