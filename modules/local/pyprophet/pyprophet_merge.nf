process PYPROPHET_MERGE {
    label 'pyprophet'
    tag   "merge"

    container "${params.pyprophet_container}"

    input:
    path osws
    path template

    output:
    path "merged.osw",       emit: merged
    path "merge.log",        emit: log
    path "versions.yml",     emit: versions

    script:
    """
    set -o pipefail

    export HOME=\$PWD
    export MPLCONFIGDIR=\$PWD/.mplconfig
    export XDG_CACHE_HOME=\$PWD/.cache
    mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME"

    pyprophet merge osw \\
        --template=${template} \\
        --out=merged.osw \\
        ${osws} \\
        2>&1 | tee merge.log

    cat <<-EOF > versions.yml
    "${task.process}":
      pyprophet: \$(pyprophet --version 2>&1 | head -n1)
    EOF
    """

    stub:
    """
    touch merged.osw merge.log
    echo '"${task.process}": {pyprophet: stub}' > versions.yml
    """
}
