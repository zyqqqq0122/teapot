process PYPROPHET_APPLY {
    label 'pyprophet'
    tag   "${meta.id}"

    container "${params.pyprophet_container}"

    publishDir path: { "${params.outdir}/pyprophet/scored/${meta.id}" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(reference_list), path(osw)
    path  model

    output:
    tuple val(meta), path(reference_list), path("${meta.id}.scored.osw"), emit: scored
    tuple val(meta), path("${meta.id}.apply.log"),                    emit: log
    path  "versions.yml", emit: versions

    script:
    def lvl = params.pyprophet_level ?: 'ms2'
    """
    set -o pipefail
    
    export HOME=\$PWD
    export MPLCONFIGDIR=\$PWD/.mplconfig
    export XDG_CACHE_HOME=\$PWD/.cache
    mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME"

    cp ${osw} ${meta.id}.scored.osw
    pyprophet score \\
        --in=${meta.id}.scored.osw \\
        --level=${lvl} \\
        --apply_weights=${model} \\
        ${Args.flat(params.pyprophet_apply_args)} \\
        2>&1 | tee ${meta.id}.apply.log

    cat <<-EOF > versions.yml
    "${task.process}":
      pyprophet: \$(pyprophet --version 2>&1 | head -n1)
    EOF
    """

    stub:
    """
    touch ${meta.id}.scored.osw ${meta.id}.apply.log
    echo '"${task.process}": {pyprophet: stub}' > versions.yml
    """
}
