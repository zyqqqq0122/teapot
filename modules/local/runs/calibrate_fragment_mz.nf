process CALIBRATE_FRAGMENT_MZ {
    label 'python'
    tag   "${meta.id}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/qc/${meta.id}" },
               mode: params.publish_mode,
               saveAs: { fn -> fn ==~ /.*\.fragment_mz\.(json|txt|targets\.tsv|calibrants\.tsv)/ ? fn : null }

    input:
    tuple val(meta), path(run), path(reference_list), path(windows, stageAs: 'windows/*'), val(fragment)
    path  library, stageAs: 'library/*'
    val   correct

    output:
    tuple val(meta), path("calibrated/${run.name}"), path(reference_list), emit: run
    tuple val(meta), path("${meta.id}.fragment_mz.targets.tsv"), optional: true, emit: located
    tuple val(meta), path("${meta.id}.fragment_mz.tolerance.tsv"), optional: true, emit: tolerance
    path  "${meta.id}.fragment_mz.row.tsv",                                emit: row
    path  "${meta.id}.fragment_mz.warnings.txt",                           emit: warnings
    path  "${meta.id}.fragment_mz.*",                                      emit: report
    path  "versions.yml",                                                   emit: versions

    script:
    def ref = reference_list.name != 'NO_FILE' ? "--reference-list ${reference_list}" : ''
    def win = windows.name.endsWith('NO_FILE') ? '' : "--windows ${windows}"
    def tol = fragment ? "--fragment-tolerance ${fragment[0]} --fragment-unit ${fragment[1]}" : ''
    """
    # bin/calibrate_fragment_mz.py ${Fingerprint.of("${projectDir}/bin/calibrate_fragment_mz.py")}: a changed script reruns the step
    calibrate_fragment_mz.py \\
        --run        ${run} \\
        --library    ${library} \\
        ${ref} \\
        ${win} \\
        ${tol} \\
        --min-offset ${params.fragment_mz_calibration_min_offset} \\
        ${correct ? '' : '--measure-only'} \\
        --sample-id  ${meta.id} \\
        --report     ${meta.id}.fragment_mz \\
        --summary    ${meta.id}.fragment_mz.row.tsv \\
        --warnings   ${meta.id}.fragment_mz.warnings.txt \\
        --tolerance-counts ${meta.id}.fragment_mz.tolerance.tsv \\
        --out        calibrated/${run.name}

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
      numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
    EOF
    """

    stub:
    """
    mkdir calibrated
    ln -s "\$(readlink -f ${run})" calibrated/${run.name}
    printf 'sample_id\\tstatus\\tcorrected\\n${meta.id}\\tstub\\tno\\n' > ${meta.id}.fragment_mz.row.tsv
    touch ${meta.id}.fragment_mz.warnings.txt ${meta.id}.fragment_mz.json
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
