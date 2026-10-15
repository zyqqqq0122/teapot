process RESTRICT_BACKGROUND {
    label 'python'
    tag   "${dlib.baseName}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/prepare_library" },
               mode: params.publish_mode,
               saveAs: { fn -> fn == 'versions.yml' ? 'restrict_background.versions.yml' : fn }

    input:
    path dlib
    path peptides, stageAs: 'peptides/*'
    val  top

    output:
    path "${dlib.baseName}.restricted.dlib", emit: dlib
    path "restrict_background.log",          emit: log
    path "versions.yml",                     emit: versions

    script:
    """
    set -o pipefail
    restrict_background.py \\
        --peptides ${peptides} \\
        --dlib     ${dlib} \\
        --top      ${top} \\
        --out      ${dlib.baseName}.restricted.dlib \\
        2>&1 | tee restrict_background.log

    cat <<-EOF > versions.yml
    "${task.process}":
      python:   \$(python3 --version 2>&1 | awk '{print \$2}')
      peptides: '${peptides.name}'
      top:      ${top}
    EOF
    """

    stub:
    """
    touch ${dlib.baseName}.restricted.dlib restrict_background.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
