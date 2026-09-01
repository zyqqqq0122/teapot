process CONTEXT_MS_RUN {
    label 'context'
    tag   "${meta.id}"

    container "${params.context_container}"

    publishDir path: { "${params.outdir}/context_ms_run/${meta.id}" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(reference_tsv), path(background_tsv)

    output:
    tuple val(meta), path("${meta.id}.psm.reference.txt"),     emit: psm
    tuple val(meta), path("${meta.id}.peptide.reference.txt"), emit: peptide
    tuple val(meta), path("${meta.id}.rescored_features.tsv"), emit: rescored
    tuple val(meta), path("weights/${meta.id}.weights.txt"),   emit: weights
    tuple val(meta), path("${meta.id}.context.log"),           emit: log
    path  "versions.yml", emit: versions

    script:
    """
    set -o pipefail

    mkdir -p weights
    python3 -m context run \\
        --background ${background_tsv} \\
        --reference  ${reference_tsv} \\
        --prefix     ${meta.id} \\
        --outdir     . \\
        --engine     ${params.context_engine} \\
        --input-profile     ${params.context_input_profile} \\
        --seed-coefficients ${params.context_seed_coefficients} \\
        ${Args.flat(params.context_args)} \\
        2>&1 | tee ${meta.id}.context.log

    cat <<-EOF > versions.yml
    "${task.process}":
      context: \$(python3 -m context --version 2>&1 | tail -n1 || echo unknown)
      engine:  ${params.context_engine}
    EOF
    """

    stub:
    """
    mkdir -p weights
    touch ${meta.id}.psm.reference.txt \\
          ${meta.id}.peptide.reference.txt \\
          ${meta.id}.rescored_features.tsv \\
          weights/${meta.id}.weights.txt \\
          ${meta.id}.context.log
    echo '"${task.process}": {context: stub}' > versions.yml
    """
}
