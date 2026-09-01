process PYPROPHET_LEARN {
    label 'pyprophet'
    tag   "learn"

    container "${params.pyprophet_container}"

    publishDir path: { "${params.outdir}/pyprophet/model" },
               mode: params.publish_mode

    input:
    path merged

    output:
    path "model.osw",                           emit: model
    path "*.pdf",                 optional: true, emit: reports
    path "learn.log",                           emit: log
    path "versions.yml",                        emit: versions

    script:
    def lvl = params.pyprophet_level ?: 'ms2'
    """
    set -o pipefail

    export HOME=\$PWD
    export MPLCONFIGDIR=\$PWD/.mplconfig
    export XDG_CACHE_HOME=\$PWD/.cache
    mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME"

    cp ${merged} model.osw
    pyprophet score --in=model.osw --level=${lvl} --threads=${task.cpus} \\
        ${Args.flat(params.pyprophet_learn_args)} \\
        2>&1 | tee learn.log

    cat <<-EOF > versions.yml
    "${task.process}":
      pyprophet: \$(pyprophet --version 2>&1 | head -n1)
    EOF
    """

    stub:
    """
    touch model.osw learn.log model_ms1ms2_report.pdf
    echo '"${task.process}": {pyprophet: stub}' > versions.yml
    """
}
