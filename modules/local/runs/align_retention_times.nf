process ALIGN_RETENTION_TIMES {
    label 'python'
    tag   "${route}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/qc" },
               mode: params.publish_mode,
               saveAs: { fn -> fn.endsWith('.apex_windows.tsv') ? "rt_alignment/${route}/${file(fn).name}" :
                               fn == "rt_alignment.${route}.tsv" ? fn : null }

    input:
    path  samples
    path  located,        stageAs: 'located/*'
    path  neighbours,     stageAs: 'neighbours/*'
    path  library,        stageAs: 'library/*'
    path  reference_list, stageAs: 'ref/*'
    val   heavy_label
    val   route

    output:
    path "windows/*.apex_windows.tsv",      emit: windows
    path "rt_alignment.${route}.tsv",       emit: summary
    path "rt_alignment.warnings.txt",       emit: warnings
    path "versions.yml",                    emit: versions

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
    # bin/align_retention_times.py ${Fingerprint.of("${projectDir}/bin/align_retention_times.py")}: a changed script reruns the step
    align_retention_times.py \\
        --samples     ${samples} \\
        --located     located \\
        --neighbours  neighbours/* \\
        --library     library/* \\
        --reference   ref/* \\
        --labels      '${label_json}' \\
        --half-window ${params.rt_alignment_half_window_seconds} \\
        --out-dir     windows \\
        --summary     rt_alignment.${route}.tsv \\
        --warnings    rt_alignment.warnings.txt

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
      numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
    EOF
    """

    stub:
    """
    mkdir windows
    for s in \$(cat ${samples}); do
        printf 'PeptideModSeq\\tPrecursorCharge\\tMinApexRTInSeconds\\tMaxApexRTInSeconds\\n' > windows/\$s.apex_windows.tsv
    done
    printf 'sample_id\\tstatus\\n' > rt_alignment.${route}.tsv
    touch rt_alignment.warnings.txt
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
