process SEARCH_GPF {
    label 'encyclopedia'
    tag   "${meta.condition}"

    publishDir path: { "${params.outdir}/search_gpf/${meta.condition}" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(combined), path(fasta)
    path  dlib

    output:
    tuple val(meta), path("${combined.name}*", includeInputs: true), emit: products
    tuple val(meta), path("${meta.condition}.search_gpf.log"), emit: log
    path  "versions.yml", emit: versions

    script:
    """
    java -Xmx${params.java_mem} -jar ${params.encyclopedia_jar} \\
        -i ${combined} \\
        -l ${dlib} \\
        -f ${fasta} \\
        ${params.gpf_search_args} \\
        2>&1 | tee ${meta.condition}.search_gpf.log

    cat <<-EOF > versions.yml
    "${task.process}":
      encyclopedia:      \$(java -jar ${params.encyclopedia_jar} -version 2>&1 | head -n1)
      condition:         ${meta.condition}
      gpf_search_args:   '${params.gpf_search_args}'
    EOF
    """

    stub:
    """
    touch ${combined.name}.features.txt
    touch ${combined.name}.encyclopedia.txt
    touch ${combined.name}.elib
    touch ${meta.condition}.search_gpf.log
    echo '"${task.process}": {encyclopedia: stub}' > versions.yml
    """
}
