process ASSERT_OPENSWATH_FEASIBLE {
    label 'python'
    tag   "${meta.id}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/openswath_feasibility" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(mzml)

    output:
    tuple val(meta), path("${meta.id}.openswath_feasibility.txt"),  emit: report
    path "versions.yml",                                            emit: versions

    script:
    def warn_only = params.openswath_feasibility_warn_only ? '--warn-only' : ''
    def width     = params.openswath_occupancy_window_width
                        ? "--window-width ${params.openswath_occupancy_window_width}" : ''
    """
    assert_openswath_feasible.py \\
        --mzml ${mzml} \\
        --tolerance ${params.openswath_occupancy_tolerance} \\
        --max-occupancy ${params.openswath_max_occupancy} \\
        --report ${meta.id}.openswath_feasibility.txt \\
        ${width} ${warn_only}

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
    EOF
    """

    stub:
    """
    touch ${meta.id}.openswath_feasibility.txt
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
