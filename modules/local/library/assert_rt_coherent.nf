process ASSERT_RT_COHERENT {
    label 'python'
    tag   "${union_dlib.baseName}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/merge_libraries" },
               mode: params.publish_mode

    input:
    path union_dlib

    output:
    path "rt_coherence.tsv", emit: report
    path "rt_coherence.log", emit: log
    path "versions.yml",     emit: versions

    script:
    def require_iqr = params.rt_coherence_require_iqr_overlap ? 1 : 0
    """
    set -o pipefail

    assert_rt_coherent.py \\
        --union ${union_dlib} \\
        --report rt_coherence.tsv \\
        --max-range-ratio ${params.rt_coherence_max_range_ratio} \\
        --require-iqr-overlap ${require_iqr} \\
        2>&1 | tee rt_coherence.log

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
    EOF
    """

    stub:
    """
    touch rt_coherence.tsv rt_coherence.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
