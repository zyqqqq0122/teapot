process QUANTIFY_STANDARD_PEAK {
    label 'python'
    tag   "${meta.id}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/quant/${route}/standard_peak" },
               mode: params.publish_mode,
               saveAs: { fn -> fn.endsWith('.standard_peak.tsv') ? fn : null }

    input:
    tuple val(meta), path(run), path(identifications, stageAs: 'ids/*')
    path  library,   stageAs: 'library/*'
    path  pairs,     stageAs: 'pairs/*'
    val   heavy_label
    val   search_args
    val   route

    output:
    path "${meta.id}.standard_peak.tsv", emit: rows
    path "versions.yml",                 emit: versions

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
    # bin/quantify_standard_peak.py ${Fingerprint.of("${projectDir}/bin/quantify_standard_peak.py")}: a changed script reruns the step
    quantify_standard_peak.py \\
        --run         ${run} \\
        --sample-id   ${meta.id} \\
        --library     ${library} \\
        --reference   ${pairs} \\
        --labels      '${label_json}' \\
        --search-args '${search_args ?: ''}' \\
        --search      ${identifications} \\
        --fdr         ${params.fdr} \\
        --out         ${meta.id}.standard_peak.tsv

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
      numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
    EOF
    """

    stub:
    """
    printf 'sample_id\\tpeptide\\tpeptidoform\\tcharge\\tstatus\\tstart\\tstop\\tlight_area\\theavy_area\\n' > ${meta.id}.standard_peak.tsv
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
