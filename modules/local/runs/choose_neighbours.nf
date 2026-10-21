process CHOOSE_NEIGHBOURS {
    label 'python'
    tag   "${route}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/qc" },
               mode: params.publish_mode,
               saveAs: { fn -> fn == "neighbours.${route}.tsv" ? fn : null }

    input:
    path  samples
    path  located, stageAs: 'located/*'
    path  library, stageAs: 'library/*'
    val   heavy_label
    val   route

    output:
    path "neighbours.${route}.tsv", emit: neighbours
    path "versions.yml",            emit: versions

    script:
    def _rows_uni = file("${projectDir}/assets/unimod_delta_masses.tsv").readLines()
    def _delta_by_acc = [:]
    _rows_uni.drop(1).each { l ->
        def c = l.split('\t')
        if (c.size() >= 3) _delta_by_acc[c[0]] = c[2] as Double
    }
    def _aug = (heavy_label ?: []).collect { r ->
        def d = r.containsKey('delta') ? (r.delta as Double) : _delta_by_acc[r.unimod]
        [ site: r.site, unimod: r.unimod, delta: d ]
    }
    def label_json = groovy.json.JsonOutput.toJson(_aug)
    """
    # bin/choose_neighbours.py ${Fingerprint.of("${projectDir}/bin/choose_neighbours.py")}: a changed script reruns the step
    # bin/align_retention_times.py ${Fingerprint.of("${projectDir}/bin/align_retention_times.py")}: a changed script reruns the step
    choose_neighbours.py \\
        --samples ${samples} \\
        --located located \\
        --library library/* \\
        --labels  '${label_json}' \\
        --out     neighbours.${route}.tsv

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
      numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
    EOF
    """

    stub:
    """
    printf 'sample_id\\tanchors\\tneighbours\\tagreement_seconds\\tpredicted_seconds\\tstatus\\n' > neighbours.${route}.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
