process CHECK_REFERENCE_COVERAGE {
    label 'python'
    tag   "${route}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/qc" },
               mode: params.publish_mode,
               saveAs: { fn -> fn.endsWith('.reference_coverage.tsv') ? fn : null }

    input:
    path reference_list
    path library
    val  route

    output:
    path "${route}.reference_coverage.tsv",          emit: coverage
    path "${route}.reference_coverage.warnings.txt", emit: warnings
    path "versions.yml",                             emit: versions

    script:
    """
    # bin/check_reference_coverage.py ${Fingerprint.of("${projectDir}/bin/check_reference_coverage.py")}: a changed script reruns the step
    set -o pipefail

    check_reference_coverage.py \\
        --reference ${reference_list} \\
        --library   ${library} \\
        --out       ${route}.reference_coverage.tsv \\
        --warnings  ${route}.reference_coverage.warnings.txt \\
        2>&1 | tee ${route}.reference_coverage.log

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
    EOF
    """

    stub:
    """
    touch ${route}.reference_coverage.tsv ${route}.reference_coverage.warnings.txt
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
