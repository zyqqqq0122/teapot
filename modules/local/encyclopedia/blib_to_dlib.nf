process BLIB_TO_DLIB {
    label 'encyclopedia'
    tag   "${blib.baseName}"

    input:
    path blib
    path fasta

    output:
    path "${blib.baseName}.dlib",              emit: dlib
    path "${blib.baseName}.blib_to_dlib.log",  emit: log
    path "versions.yml",                       emit: versions

    script:
    // -convert -blibToLib (ConvertBLIBToLibrary.java)
    //   -i <blib> -f <fasta>  required
    //   -o <dlib>             output (optional; defaults to <blib>.dlib)
    """
    java -Xmx${params.java_mem} -jar ${params.encyclopedia_jar} \\
        -convert -blibToLib \\
        -i ${blib} \\
        -f ${fasta} \\
        -o ${blib.baseName}.dlib \\
        2>&1 | tee ${blib.baseName}.blib_to_dlib.log

    cat <<-EOF > versions.yml
    "${task.process}":
      encyclopedia: \$(java -jar ${params.encyclopedia_jar} -version 2>&1 | head -n1)
    EOF
    """

    stub:
    """
    touch ${blib.baseName}.dlib
    touch ${blib.baseName}.blib_to_dlib.log
    echo '"${task.process}": {encyclopedia: stub}' > versions.yml
    """
}
