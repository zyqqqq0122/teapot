process BUILD_CHRLIB {
    label 'encyclopedia'
    tag   "elib"

    input:
    path search_products
    path dlib
    path fasta

    output:
    path "chromatogram.elib", emit: elib
    path "build_chrlib.log",  emit: log
    path "versions.yml",      emit: versions

    script:
    // -libexport aggregation.
    //   -i .   directory scan; picks up every .dia + features.txt + .elib in cwd
    //   -l     original searched library (.dlib or .elib)
    //   -f     fasta used across the searches
    //   ${params.libexport_args}  pass-through knobs (e.g. -a false, -blib, ...)
    """
    exec > >(tee -a build_chrlib.log) 2>&1

    echo "==> libexport aggregation on \$(ls -1 | wc -l) staged files"

    java -Xmx${params.java_mem} -jar ${params.encyclopedia_jar} \\
        -libexport \\
        -o chromatogram.elib \\
        -i . \\
        -l ${dlib} \\
        -f ${fasta} \\
        ${params.libexport_args}

    cat <<-EOF > versions.yml
    "${task.process}":
      encyclopedia:    \$(java -jar ${params.encyclopedia_jar} -version 2>&1 | head -n1)
      libexport_args:  '${params.libexport_args}'
    EOF
    """

    stub:
    """
    touch chromatogram.elib build_chrlib.log
    echo '"${task.process}": {encyclopedia: stub}' > versions.yml
    """
}
