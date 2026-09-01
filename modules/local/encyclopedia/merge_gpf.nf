process MERGE_GPF {
    label 'encyclopedia'
    tag   "${meta.condition}"

    publishDir path: { "${params.outdir}/merge_gpf/${meta.condition}" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(fractions)

    output:
    tuple val(meta), path("${meta.condition}_combined.dia"),   emit: combined
    tuple val(meta), path("${meta.condition}.merge_gpf.log"),  emit: log
    path  "versions.yml", emit: versions

    script:
    def joined = fractions.collect { it.name }.join(':')
    """
    set -o pipefail

    java -Xmx${params.java_mem} -jar ${params.encyclopedia_jar} \\
        -convert -mergeDIA \\
        -i ${joined} \\
        -o ${meta.condition}_combined.dia \\
        2>&1 | tee ${meta.condition}.merge_gpf.log

    cat <<-EOF > versions.yml
    "${task.process}":
      encyclopedia: \$(java -jar ${params.encyclopedia_jar} -version 2>&1 | head -n1)
      condition:    ${meta.condition}
      n_fractions:  ${fractions.size()}
    EOF
    """

    stub:
    """
    touch ${meta.condition}_combined.dia ${meta.condition}.merge_gpf.log
    echo '"${task.process}": {encyclopedia: stub}' > versions.yml
    """
}
