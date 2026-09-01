process PYPROPHET_EXPORT {
    label 'pyprophet'
    tag   "${meta.id}"

    container "${params.pyprophet_container}"

    publishDir path: { "${params.outdir}/pyprophet/export/${meta.id}" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(reference_list), path(osw)

    output:
    tuple val(meta), path(reference_list), path("${meta.id}.tsv"),       emit: tsv
    tuple val(meta), path("${meta.id}.export.log"),                 emit: log
    path  "versions.yml", emit: versions

    script:
    """
    set -o pipefail
    
    export HOME=\$PWD
    export MPLCONFIGDIR=\$PWD/.mplconfig
    export XDG_CACHE_HOME=\$PWD/.cache
    mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME"

    pyprophet export tsv \\
        --in=${osw} \\
        --out=${meta.id}.tsv \\
        --format=legacy_merged \\
        --max_rs_peakgroup_qvalue=${params.fdr} \\
        ${Args.flat(params.pyprophet_export_args)} \\
        2>&1 | tee ${meta.id}.export.log

    cat <<-EOF > versions.yml
    "${task.process}":
      pyprophet: \$(pyprophet --version 2>&1 | head -n1)
    EOF
    """

    stub:
    """
    touch ${meta.id}.tsv ${meta.id}.export.log
    echo '"${task.process}": {pyprophet: stub}' > versions.yml
    """
}
