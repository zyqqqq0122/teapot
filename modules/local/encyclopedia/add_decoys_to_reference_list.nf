process ADD_DECOYS_TO_REFERENCE_LIST {
    label 'encyclopedia'
    tag   "${reference_list.baseName}"

    publishDir path: { "${params.outdir}/prepare_library" },
               mode: params.publish_mode

    input:
    path reference_list
    path decoy_jar

    output:
    path "${reference_list.baseName}.with_decoys.txt", emit: reference_list
    path "${reference_list.baseName}.add_decoys.log",  emit: log
    path "versions.yml",                               emit: versions

    script:
    """
    set -o pipefail

    java -Xmx${params.java_mem} \\
        -cp ${params.encyclopedia_jar}:${decoy_jar} \\
        MassListDecoyGenerator \\
        -massList ${reference_list} \\
        -o ${reference_list.baseName}.with_decoys.txt \\
        2>&1 | tee ${reference_list.baseName}.add_decoys.log

    cat <<-EOF > versions.yml
    "${task.process}":
      encyclopedia: \$(java -jar ${params.encyclopedia_jar} -version 2>&1 | head -n1)
      generator: MassListDecoyGenerator (PeptideUtils.getSmartDecoy)
    EOF
    """

    stub:
    """
    touch ${reference_list.baseName}.with_decoys.txt ${reference_list.baseName}.add_decoys.log
    echo '"${task.process}": {encyclopedia: stub}' > versions.yml
    """
}
