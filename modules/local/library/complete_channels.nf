process COMPLETE_CHANNELS {
    label 'python'
    tag   "${library.baseName}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/complete_channels" },
               mode: params.publish_mode

    input:
    path library
    val  heavy_label

    output:
    path "${library.baseName}.complete.${library.extension}", emit: library
    path "complete_channels_qc.tsv",                            emit: qc
    path "complete_channels.log",                               emit: log
    path "versions.yml",                                       emit: versions

    script:
    def _unimod = file("${projectDir}/assets/unimod_delta_masses.tsv").readLines()
    def _delta_by_acc = [:]
    _unimod.drop(1).each { l ->
        def c = l.split('\t')
        if (c.size() >= 3) _delta_by_acc[c[0]] = c[2] as Double
    }
    def _aug = (heavy_label ?: []).collect { r ->
        def d = r.containsKey('delta') ? (r.delta as Double) : _delta_by_acc[r.unimod]
        [ site: (r.site ?: r.residue), unimod: (r.unimod ?: ''), delta: d ]
    }
    def label_json = groovy.json.JsonOutput.toJson(_aug)
    def ext        = library.extension
    def out_name   = "${library.baseName}.complete.${ext}"
    def ion_types  = (params.complete_channels_ion_types as List).join(',')
    """
    set -o pipefail

    complete_channels.py \\
        --in ${library} \\
        --out ${out_name} \\
        --labels '${label_json}' \\
        --mod-tol ${params.mod_mass_tolerance} \\
        --ppm ${params.complete_channels_tolerance_ppm} \\
        --ion-types '${ion_types}' \\
        --max-frag-charge ${params.complete_channels_max_frag_charge} \\
        --coverage-hard-floor ${params.complete_channels_coverage_hard_floor} \\
        --coverage-warn-below ${params.complete_channels_coverage_warn_below} \\
        --skip-proteins '${(params.complete_channels_skip_proteins ?: []).join(',')}' \\
        --qc complete_channels_qc.tsv \\
        2>&1 | tee complete_channels.log

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
      labels_resolved: '${label_json}'
    EOF
    """

    stub:
    """
    touch ${library.baseName}.complete.${library.extension} \\
          complete_channels_qc.tsv complete_channels.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
