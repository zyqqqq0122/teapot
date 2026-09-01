process EMIT_IRT_LIBRARY {
    label 'python'
    tag   "${blib.baseName}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/emit_irt_library" },
               mode: params.publish_mode

    input:
    path blib
    val  heavy_label

    output:
    path "${blib.baseName}.irt_anchors.tsv",  emit: anchors
    path "${blib.baseName}.irt_library.tsv",  emit: osw_tsv
    path "emit_irt_library.log",              emit: log
    path "versions.yml",                      emit: versions

    script:
    def _unimod = file("${projectDir}/assets/unimod_delta_masses.tsv").readLines()
    def _delta_by_acc = [:]
    _unimod.drop(1).each { l ->
        def c = l.split('\t')
        if (c.size() >= 3) _delta_by_acc[c[0]] = c[2] as Double
    }
    def _deltas = (heavy_label ?: []).collect { r ->
        r.containsKey('delta') ? (r.delta as Double) : _delta_by_acc[r.unimod]
    }.findAll { it != null }
    def deltas_arg = _deltas ? _deltas.collect { String.format('%.4f', it) }.join(',') : ''
    """
    set -o pipefail

    emit_irt_library.py \\
        --in ${blib} \\
        --anchors-tsv ${blib.baseName}.irt_anchors.tsv \\
        --osw-tsv     ${blib.baseName}.irt_library.tsv \\
        --top-transitions ${params.irt_top_transitions} \\
        --require-label-deltas '${deltas_arg}' \\
        2>&1 | tee emit_irt_library.log

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
      required_label_deltas: '${deltas_arg}'
    EOF
    """

    stub:
    """
    touch ${blib.baseName}.irt_anchors.tsv \\
          ${blib.baseName}.irt_library.tsv \\
          emit_irt_library.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
