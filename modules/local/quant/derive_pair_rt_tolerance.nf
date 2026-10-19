process DERIVE_PAIR_RT_TOLERANCE {
    label 'python'
    tag   "${route}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/quant/${route}" },
               mode: params.publish_mode,
               saveAs: { fn -> fn == 'versions.yml' ? 'derive_pair_rt_tolerance.versions.yml' : fn }

    input:
    path base_long
    val  route

    output:
    path "pair_rt_tolerance.tsv", emit: tolerance
    path "versions.yml",          emit: versions

    script:
    """
    # bin/derive_pair_rt_tolerance.py ${Fingerprint.of("${projectDir}/bin/derive_pair_rt_tolerance.py")}: a changed script reruns the step
    derive_pair_rt_tolerance.py --base ${base_long} --primary ${params.primary_abundance} --fdr ${params.fdr} \\
        --setting '${params.pair_rt_tolerance_seconds}' --out pair_rt_tolerance.tsv

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
    EOF
    """

    stub:
    """
    printf 'setting\\ttolerance_seconds\\tsource\\n${params.pair_rt_tolerance_seconds}\\t9\\tstub\\n' > pair_rt_tolerance.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
