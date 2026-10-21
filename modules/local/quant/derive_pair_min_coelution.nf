process DERIVE_PAIR_MIN_COELUTION {
    label 'python'
    tag   "${route}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/quant/${route}" },
               mode: params.publish_mode,
               saveAs: { fn -> fn == 'versions.yml' ? 'derive_pair_min_coelution.versions.yml' : fn }

    input:
    path tables, stageAs: 'standard_peak/*'
    path neighbours, stageAs: 'neighbours/*'
    val  route

    output:
    path "pair_min_coelution.tsv", emit: min_coelution
    path "versions.yml",           emit: versions

    script:
    def nb = neighbours.name == 'NO_FILE' ? '' : "--neighbours ${neighbours}"
    """
    # bin/derive_pair_min_coelution.py ${Fingerprint.of("${projectDir}/bin/derive_pair_min_coelution.py")}: a changed script reruns the step
    derive_pair_min_coelution.py --tables standard_peak ${nb} --setting '${params.pair_min_coelution}' \\
        --out pair_min_coelution.tsv

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
      numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
    EOF
    """

    stub:
    """
    printf 'sample_id\\tsetting\\tmin_coelution\\tchance_pass\\tmoved_traces\\tpooled_with\\tsource\\n*\\t${params.pair_min_coelution}\\t0.7\\t\\t0\\t\\tstub\\n' > pair_min_coelution.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
