process ENCYCLOPEDIA_LIBEXPORT {
    label 'encyclopedia'
    tag   "${output_basename}"

    publishDir path: { "${params.outdir}/encyclopedia_libexport/${output_basename}" },
               mode: params.publish_mode

    input:
    path search_products
    val  output_basename
    path library
    path fasta
    val  search_args

    output:
    path "${output_basename}.elib",                        emit: elib
    path "${output_basename}.elib.peptides.txt", optional: true, emit: peptides
    path "${output_basename}.elib.proteins.txt", optional: true, emit: proteins
    path "${output_basename}.libexport.log",               emit: log
    path "versions.yml",                                    emit: versions

    script:
    def jar_fp = Fingerprint.of(params.encyclopedia_jar)
    def inherited = Args.inherit(search_args)
    def extra     = Args.flat(params.libexport_args)
    """
    exec > >(tee -a ${output_basename}.libexport.log) 2>&1

    echo "==> libexport (${output_basename}) on \$(ls -1 | wc -l) staged files"
    echo "==> encyclopedia jar: ${params.encyclopedia_jar} [${jar_fp}]"
    echo "==> inherited search args: ${inherited ?: '(none)'}"
    echo "==> libexport_args:        ${extra ?: '(none)'}"

    java -Xmx${params.java_mem} -cp ${params.encyclopedia_jar} \\
        edu.washington.gs.maccoss.encyclopedia.Encyclopedia \\
        -libexport \\
        -o ${output_basename}.elib \\
        -i . \\
        -l ${library} \\
        -f ${fasta} \\
        ${inherited} ${extra}

    cat <<-EOF > versions.yml
    "${task.process}":
      encyclopedia:    \$(java -jar ${params.encyclopedia_jar} -version 2>&1 | head -n1)
      basename:        '${output_basename}'
      inherited_args:  '${inherited}'
      libexport_args:  '${extra}'
    EOF
    """

    stub:
    """
    touch ${output_basename}.elib \\
          ${output_basename}.elib.peptides.txt \\
          ${output_basename}.elib.proteins.txt \\
          ${output_basename}.libexport.log
    echo '"${task.process}": {encyclopedia: stub}' > versions.yml
    """
}
