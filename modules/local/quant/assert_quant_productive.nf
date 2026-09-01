process ASSERT_QUANT_PRODUCTIVE {
    label 'python'
    tag   "${engine}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/qc" },
               mode: params.publish_mode

    input:
    path  analysis_elib
    val   engine

    output:
    path "${engine}.quant_qc.tsv", emit: qc
    path "versions.yml",           emit: versions

    script:
    """
    python3 - <<'PY' > ${engine}.quant_qc.tsv
    import sqlite3, sys, csv

    ELIB     = "${analysis_elib}"
    ENGINE   = "${engine}"
    MAX_ZERO = float("${params.quant_max_zero_ion_fraction}")
    MIN_RUNS = int("${params.quant_min_runs}")

    db = sqlite3.connect(ELIB)
    c  = db.cursor()

    have = {r[0] for r in c.execute("select name from sqlite_master where type='table'")}
    if 'peptidequants' not in have:
        sys.exit("ASSERT_QUANT_PRODUCTIVE [%s]: %s has no peptidequants table, so "
                 "libexport quantified nothing at all." % (ENGINE, ELIB))

    total, zero_ions, zero_int, mean_ions, runs = c.execute(
        "select count(*), sum(NumberOfQuantIons=0), sum(TotalIntensity=0), "
        "       avg(NumberOfQuantIons), count(distinct SourceFile) "
        "from peptidequants").fetchone()

    zero_ions = zero_ions or 0
    zero_int  = zero_int or 0
    frac      = (zero_ions / total) if total else 1.0

    w = csv.writer(sys.stdout, delimiter='\\t')
    w.writerow(['engine', 'peptidequant_rows', 'source_files', 'zero_quant_ions',
                'zero_intensity', 'zero_ion_fraction', 'mean_quant_ions'])
    w.writerow([ENGINE, total, runs, zero_ions, zero_int,
                '%.4f' % frac, '%.3f' % (mean_ions or 0.0)])
    sys.stdout.flush()

    if total == 0:
        sys.exit("ASSERT_QUANT_PRODUCTIVE [%s]: peptidequants is empty -- libexport "
                 "ran but extracted no peptide at all." % ENGINE)

    if zero_ions == total:
        sys.exit(
            "ASSERT_QUANT_PRODUCTIVE [%s]: NumberOfQuantIons is 0 for all %d "
            "peptides, so every abundance in this run is 0. The search identified "
            "these peptides; the quant pass then matched no fragment ion for any "
            "of them.\\n"
            "  The usual cause is that libexport extracted at a DIFFERENT mass "
            "tolerance from the search. libexport parses its own command line and "
            "defaults to Orbitrap (-ftol 10 -ftolunits ppm, -frag CID); on ion-trap "
            "MS2 the fragment error is ~0.1-0.3 Da, so nothing matches and every "
            "chromatogram is zero.\\n"
            "  The pipeline inherits -ftol/-frag from the search args for you -- "
            "check the 'inherited search args' line at the top of the libexport "
            "log. If it says (none), set dia_search_args / context_search_args / "
            "gpf_search_args for this run, e.g. '-ftol 0.4 -ftolunits AMU "
            "-lftol 0.4 -lftolunits AMU -frag HCD'." % (ENGINE, total))

    if frac > MAX_ZERO:
        print("WARNING: ASSERT_QUANT_PRODUCTIVE [%s]: %d of %d peptides (%.0f%%) have "
              "NumberOfQuantIons=0 and therefore zero abundance, above the %.0f%% "
              "threshold (quant_max_zero_ion_fraction). Some are legitimately absent "
              "from this run, but check the 'inherited search args' line in the "
              "libexport log before trusting the abundances."
              % (ENGINE, zero_ions, total, 100 * frac, 100 * MAX_ZERO), file=sys.stderr)

    if runs < MIN_RUNS:
        print("WARNING: ASSERT_QUANT_PRODUCTIVE [%s]: quantification is based on %d "
              "run(s), below quant_min_runs=%d. EncyclopeDIA needs at least two runs "
              "to align retention times and write a cross-run quant matrix, and its "
              "abundances are calibrated across a cohort -- with this few runs treat "
              "them as relative within-run values only."
              % (ENGINE, runs, MIN_RUNS), file=sys.stderr)
    PY

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
    EOF
    """

    stub:
    """
    echo -e "engine\\tpeptidequant_rows\\tsource_files\\tzero_quant_ions\\tzero_intensity\\tzero_ion_fraction\\tmean_quant_ions" > ${engine}.quant_qc.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
