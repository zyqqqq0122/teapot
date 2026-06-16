process PYPROPHET {
    tag "${sample}"
    label 'pyprophet'
    publishDir "${params.outdir}/pyprophet/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(osw)

    output:
    tuple val(sample), path("${sample}_scored.osw"), emit: scored

    script:
    """
    pyprophet score \\
        --in=${osw} \\
        --out=${sample}_scored.osw \\
        --level=ms2
    """
}
