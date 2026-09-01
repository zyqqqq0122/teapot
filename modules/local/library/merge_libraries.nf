process MERGE_LIBRARIES {
    label 'python'
    tag   "${primary.baseName}+${secondary.baseName}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/merge_libraries" },
               mode: params.publish_mode

    input:
    path primary
    path secondary

    output:
    path "${primary.baseName}.union.dlib", emit: dlib
    path "merge_libraries.log",              emit: log
    path "versions.yml",                   emit: versions

    script:
    """
    set -o pipefail

    merge_libraries.py \\
        --primary ${primary} \\
        --secondary ${secondary} \\
        --out ${primary.baseName}.union.dlib \\
        2>&1 | tee merge_libraries.log

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
    EOF
    """

    stub:
    """
    touch ${primary.baseName}.union.dlib merge_libraries.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
