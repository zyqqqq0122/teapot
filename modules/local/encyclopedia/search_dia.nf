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
    // features:     raw scored features across all candidates
    // targets:      Percolator-filtered target peptides at params.fdr (canonical DIA result)
    // decoys:       Percolator-filtered decoys, useful for QC / FDR sanity checks
    // elib:         per-file chromatogram-library output written by encyclopedia
    tuple val(meta), path("${dia.baseName}.features.txt"),                              emit: features
    tuple val(meta), path("${dia.baseName}.encyclopedia.txt"),       optional: true,     emit: targets
    tuple val(meta), path("${dia.baseName}.encyclopedia.decoy.txt"), optional: true,     emit: decoys
    tuple val(meta), path("${dia.baseName}.elib"),                    optional: true,     emit: elib
    tuple val(meta), path("${meta.id}.search.log"),                                       emit: log
    path  "versions.yml", emit: versions

    script:
    """
    java -Xmx${params.java_mem} -jar ${params.encyclopedia_jar} \\
        -i ${dia} \\
        -l ${library} \\
        -f ${fasta} \\
        ${params.dia_search_args} \\
        2>&1 | tee ${meta.id}.search.log

    cat <<-EOF > versions.yml
    "${task.process}":
      encyclopedia: \$(java -jar ${params.encyclopedia_jar} -version 2>&1 | head -n1)
    EOF
    """

    stub:
    """
    touch ${dia.baseName}.features.txt \\
          ${dia.baseName}.encyclopedia.txt \\
          ${dia.baseName}.encyclopedia.decoy.txt \\
          ${dia.baseName}.elib \\
          ${meta.id}.search.log
    echo '"${task.process}": {encyclopedia: stub}' > versions.yml
    """
}
