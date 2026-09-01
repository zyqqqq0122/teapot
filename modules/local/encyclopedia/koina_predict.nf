process KOINA_PREDICT {
    label 'encyclopedia'
    tag   "${fasta.baseName}"

    input:
    path fasta

    output:
    path "${fasta.baseName}.koina.dlib", emit: dlib
    path "${fasta.baseName}.koina.log",  emit: log
    path "versions.yml",                 emit: versions

    script:
    """
    set -o pipefail

    java -Xmx${params.java_mem} -jar ${params.encyclopedia_jar} \\
        -convert -fastaToKoinaLibrary \\
        -i ${fasta} \\
        -o ${fasta.baseName}.koina.dlib \\
        -url '${params.koina_url}' \\
        -models '${params.koina_models}' \\
        2>&1 | tee ${fasta.baseName}.koina.log

    cat <<-EOF > versions.yml
    "${task.process}":
      encyclopedia: \$(java -jar ${params.encyclopedia_jar} -version 2>&1 | head -n1)
      koina_models: '${params.koina_models}'
      koina_url:    '${params.koina_url}'
    EOF
    """

    stub:
    """
    touch ${fasta.baseName}.koina.dlib
    touch ${fasta.baseName}.koina.log
    echo '"${task.process}": {encyclopedia: stub, koina_models: "${params.koina_models}"}' > versions.yml
    """
}
