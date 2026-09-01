process ENCYCLOPEDIA_LIBEXPORT_PER_SAMPLE {
    label 'encyclopedia'
    tag   "${meta.id}"

    publishDir path: { "${params.outdir}/encyclopedia_libexport_per_sample/${meta.id}" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(products)
    path  library
    path  fasta
    val   search_args

    output:
    tuple val(meta), path("*.elib", includeInputs: false), emit: elib
    tuple val(meta), path("${meta.id}.libexport_single.log"), emit: log
    path  "versions.yml", emit: versions

    script:
    def jar_fp    = Fingerprint.of(params.encyclopedia_jar)
    def inherited = Args.inherit(search_args)
    def extra     = Args.flat(params.libexport_args)
    """
    exec > >(tee -a ${meta.id}.libexport_single.log) 2>&1

    echo "==> per-sample libexport for ${meta.id}"
    echo "==> inherited search args: ${inherited ?: '(none)'}"

    src=\$(ls *.encyclopedia.txt | head -n1)
    stem=\${src%.encyclopedia.txt}
    echo "==> source: \$src  ->  \${stem}.elib"

    java -Xmx${params.java_mem} -cp ${params.encyclopedia_jar} \\
        edu.washington.gs.maccoss.encyclopedia.Encyclopedia \\
        -libexport \\
        -o "\${stem}.elib" \\
        -i . \\
        -l ${library} \\
        -f ${fasta} \\
        ${inherited} ${extra}

    cat <<-EOF > versions.yml
    "${task.process}":
      encyclopedia:   \$(java -jar ${params.encyclopedia_jar} -version 2>&1 | head -n1)
      inherited_args: '${inherited}'
    EOF
    """

    stub:
    """
    touch ${meta.id}.stub.elib ${meta.id}.libexport_single.log
    echo '"${task.process}": {encyclopedia: stub}' > versions.yml
    """
}
