process SEARCH_DIA {
    label 'encyclopedia'
    tag   "${meta.id}"

    publishDir path: { "${params.outdir}/search_dia/${meta.id}" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(dia)
    path  library
    path  fasta

    output:
    tuple val(meta), path("${dia.name}.features.txt"),                              emit: features
    tuple val(meta), path("${dia.name}.encyclopedia.txt"),       optional: true,     emit: targets
    tuple val(meta), path("${dia.name}.encyclopedia.decoy.txt"), optional: true,     emit: decoys
    tuple val(meta), path("${dia.name}.elib"),                    optional: true,     emit: elib
    tuple val(meta), path("${meta.id}.search.log"),                                       emit: log
    path  "versions.yml", emit: versions

    script:
    def jar_fp = Fingerprint.of(params.encyclopedia_jar)
    """
    set -o pipefail

    java -Xmx${params.java_mem} -cp ${params.encyclopedia_jar} \\
        edu.washington.gs.maccoss.encyclopedia.Encyclopedia \\
        -i ${dia} \\
        -l ${library} \\
        -f ${fasta} \\
        ${Args.flat(params.dia_search_args)} \\
        2>&1 | tee ${meta.id}.search.log

    cat <<-EOF > versions.yml
    "${task.process}":
      encyclopedia: \$(java -jar ${params.encyclopedia_jar} -version 2>&1 | head -n1)
    EOF
    """

    stub:
    """
    touch ${dia.name}.features.txt \\
          ${dia.name}.encyclopedia.txt \\
          ${dia.name}.encyclopedia.decoy.txt \\
          ${dia.name}.elib \\
          ${meta.id}.search.log
    echo '"${task.process}": {encyclopedia: stub}' > versions.yml
    """
}
