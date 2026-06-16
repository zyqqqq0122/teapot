process ENCYCLOPEDIA {
    tag "${sample}"
    label 'encyclopedia'
    publishDir "${params.outdir}/encyclopedia/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(mzml)
    path  spectral_lib
    path  fasta

    output:
    tuple val(sample), path("${sample}.elib"),  emit: elib
    tuple val(sample), path("${sample}.encyclopedia.txt"), emit: quant

    script:
    """
    encyclopedia \\
        -i ${mzml} \\
        -f ${fasta} \\
        -l ${spectral_lib} \\
        -o ${sample}.elib \\
        -numberOfThreadsUsed ${task.cpus} \\
        ${params.encyclopedia_params}
    """
}
