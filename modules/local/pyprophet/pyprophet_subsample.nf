process PYPROPHET_SUBSAMPLE {
    label 'pyprophet'
    tag   "${meta.id}"

    container "${params.pyprophet_container}"

    input:
    tuple val(meta), path(reference_list), path(osw)
    val   n_runs

    output:
    tuple val(meta), path("${meta.id}.osws"),  emit: osws
    tuple val(meta), path("${meta.id}.subsample.log"), emit: log
    path  "versions.yml", emit: versions

    script:
    def ratio = n_runs > 0 ? Math.min(1.0d, 1.0d / (n_runs as double)) : 1.0d
    """
    set -o pipefail

    export HOME=\$PWD
    export MPLCONFIGDIR=\$PWD/.mplconfig
    export XDG_CACHE_HOME=\$PWD/.cache
    mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME"

    pyprophet subsample \\
        --in=${osw} \\
        --out=${meta.id}.osws \\
        --subsample_ratio=${ratio} \\
        2>&1 | tee ${meta.id}.subsample.log

    cat <<-EOF > versions.yml
    "${task.process}":
      pyprophet: \$(pyprophet --version 2>&1 | head -n1)
    EOF
    """

    stub:
    """
    touch ${meta.id}.osws ${meta.id}.subsample.log
    echo '"${task.process}": {pyprophet: stub}' > versions.yml
    """
}
