process PYPROPHET_SCORE_SINGLE {
    label 'pyprophet'
    tag   "${meta.id}"

    container "${params.pyprophet_container}"

    publishDir path: { "${params.outdir}/pyprophet_prm/${meta.id}" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(osw)

    output:
    tuple val(meta), path("${meta.id}.scored.osw"), emit: scored
    tuple val(meta), path("${meta.id}.tsv"),        emit: tsv
    tuple val(meta), path("${meta.id}.pyprophet.log"), emit: log
    path  "*.pdf", optional: true,                  emit: reports
    path  "versions.yml",                           emit: versions

    script:
    def lvl = params.pyprophet_level ?: 'ms2'
    """
    export HOME=\$PWD
    export MPLCONFIGDIR=\$PWD/.mplconfig
    export XDG_CACHE_HOME=\$PWD/.cache
    mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME"
    set -o pipefail

    cp ${osw} ${meta.id}.scored.osw
    pyprophet score \\
        --in=${meta.id}.scored.osw \\
        --level=${lvl} \\
        --threads=${task.cpus} \\
        ${Args.flat(params.pyprophet_learn_args)} \\
        2>&1 | tee ${meta.id}.pyprophet.log

    pyprophet export tsv \\
        --in=${meta.id}.scored.osw \\
        --out=${meta.id}.tsv \\
        --format=legacy_merged \\
        --max_rs_peakgroup_qvalue=${params.fdr} \\
        ${Args.flat(params.pyprophet_export_args)} \\
        2>&1 | tee -a ${meta.id}.pyprophet.log

    cat <<-EOF > versions.yml
    "${task.process}":
      pyprophet: \$(pyprophet --version 2>&1 | head -n1)
      level: '${lvl}'
    EOF
    """

    stub:
    """
    touch ${meta.id}.scored.osw ${meta.id}.tsv ${meta.id}.pyprophet.log
    echo '"${task.process}": {pyprophet: stub}' > versions.yml
    """
}
