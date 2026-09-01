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
    tuple val(meta), path("nested_windows.flag"), optional: true,   emit: nested
    path "versions.yml",                                            emit: versions

    script:
    def warn_only = params.openswath_feasibility_warn_only ? '--warn-only' : ''
    def pinned = ("${params.openswath_search_args}" =~ /-matching_window_only/)
    def mwo    = pinned ? '--matching-window-only-set' : ''
    def auto   = (!pinned && params.openswath_auto_matching_window)
                     ? '--auto-matching-window --emit-nested-flag nested_windows.flag' : ''
    def width     = params.openswath_occupancy_window_width
                        ? "--window-width ${params.openswath_occupancy_window_width}" : ''
    """
    assert_openswath_feasible.py \\
        --mzml ${mzml} \\
        --tolerance ${params.openswath_occupancy_tolerance} \\
        --max-occupancy ${params.openswath_max_occupancy} \\
        --report ${meta.id}.openswath_feasibility.txt \\
        ${width} ${warn_only} ${mwo} ${auto}

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
