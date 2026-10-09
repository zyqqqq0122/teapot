process SELECT_DIATHEM_TARGETS {
    label 'python'
    tag   "${id_format}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/diathem_quant/DIA" },
               mode: params.publish_mode,
               saveAs: { fn -> fn.startsWith('prior/') ? null : fn }

    input:
    path ids,     stageAs: 'ids/*'
    val  id_format
    path library, stageAs: 'library/*'
    path prior,   stageAs: 'prior_in/*'
    path unimod

    output:
    path "diathem_targets.csv",        emit: targets
    path "prior/*", optional: true,    emit: prior
    path "select_diathem_targets.log", emit: log
    path "versions.yml",               emit: versions

    script:
    def lib_arg = library.name.tokenize('/').last() == 'NO_FILE' ? '' : "--library ${library}"
    """
    set -o pipefail

    select_diathem_targets.py \\
        --ids ids \\
        --format ${id_format} \\
        --fdr ${params.fdr} \\
        ${lib_arg} \\
        --prior ${prior} \\
        --unimod ${unimod} \\
        --out diathem_targets.csv \\
        --prior-out-dir prior \\
        2>&1 | tee select_diathem_targets.log

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
      fdr:    ${params.fdr}
    EOF
    """

    stub:
    """
    touch diathem_targets.csv select_diathem_targets.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
