process REMOVE_SURVEY_SCANS {
    label 'python'
    tag   "${meta.id}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/isolation_windows" },
               mode: params.publish_mode,
               pattern: '*.isolation_windows.tsv'

    input:
    tuple val(meta), path(run), path(reference_list)

    output:
    tuple val(meta), path("targeted/${run.name}"), path(reference_list), emit: run
    tuple val(meta), path("${meta.id}.isolation_windows.tsv"),           emit: windows
    path  "versions.yml",                                                 emit: versions

    script:
    """
    # bin/remove_survey_scans.py ${Fingerprint.of("${projectDir}/bin/remove_survey_scans.py")}: a changed script reruns the step
    remove_survey_scans.py \\
        --run   ${run} \\
        --out   targeted/${run.name} \\
        --table ${meta.id}.isolation_windows.tsv

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
    EOF
    """

    stub:
    """
    mkdir targeted
    ln -s "\$(readlink -f ${run})" targeted/${run.name}
    touch ${meta.id}.isolation_windows.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
