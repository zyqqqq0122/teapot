process DERIVE_SEARCH_TOLERANCE {
    label 'python'
    tag   "${where}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/qc" },
               mode: params.publish_mode,
               saveAs: { fn -> fn == 'versions.yml' ? null : fn }

    input:
    path counts, stageAs: 'counts/*'
    val  fragment
    val  where

    output:
    path "fragment_tolerance.${where}.tsv",      emit: tolerance
    path "fragment_tolerance.${where}.runs.tsv", emit: runs
    path "versions.yml",                          emit: versions

    script:
    """
    # bin/derive_search_tolerance.py ${Fingerprint.of("${projectDir}/bin/derive_search_tolerance.py")}: a changed script reruns the step
    derive_search_tolerance.py --counts counts \\
        --fragment-tolerance ${fragment[0]} --fragment-unit ${fragment[1]} \\
        --setting ${Args.mode(params.search_tolerance)} --where ${where} \\
        --out fragment_tolerance.${where}.tsv --runs-out fragment_tolerance.${where}.runs.tsv

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
      numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
    EOF
    """

    stub:
    """
    printf 'route\\tsetting\\tsearch_tolerance\\tunit\\truns_at_it\\tsource\\n${where}\\t${Args.mode(params.search_tolerance)}\\t${fragment[0]}\\t${fragment[1] == 'th' ? 'Th' : 'ppm'}\\t0\\tstub\\n' > fragment_tolerance.${where}.tsv
    printf 'sample_id\\tfragment_errors\\town_search_tolerance\\tsearched_at\\n' > fragment_tolerance.${where}.runs.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
