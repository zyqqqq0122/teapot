process OPENSWATH {
    tag "${sample}"
    label 'openswath'
    publishDir "${params.outdir}/openswath/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(mzml)
    path  spectral_lib

    output:
    tuple val(sample), path("${sample}.osw"), emit: osw
    tuple val(sample), path("${sample}_chrom.mzML"), emit: chrom

    script:
    """
    OpenSwathWorkflow \\
        -in ${mzml} \\
        -tr ${spectral_lib} \\
        -out_osw ${sample}.osw \\
        -out_chrom ${sample}_chrom.mzML \\
        -threads ${task.cpus} \\
        ${params.openswath_params}
    """
}
