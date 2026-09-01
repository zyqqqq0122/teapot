process HARMONIZE_RT {
    label 'python'
    tag   "${library.baseName}/${provenance}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/harmonize_rt/${provenance}" },
               mode: params.publish_mode

    input:
    path library
    path anchors
    val  provenance
    val  passthrough
    val  emit_full_window

    output:
    path "${library.baseName}.irt.${library.extension}",  emit: library
    path "${library.baseName}.rt_fit.tsv",                emit: fit_report
    path "harmonize_rt.log",                              emit: log
    path "rt_extraction_full_window",                     emit: full_window, optional: true
    path "versions.yml",                                  emit: versions

    script:
    def out_name    = "${library.baseName}.irt.${library.extension}"
    def report_name = "${library.baseName}.rt_fit.tsv"
    def pass_flag   = passthrough ? '--assume-koina-passthrough' : ''
    def anchors_arg = passthrough ? '' : "--anchors ${anchors}"
    def fallback    = params.rt_harmonize_fallback ?: 'fail'
    def full_wnd    = emit_full_window ? '--emit-full-window-sentinel' : ''
    def cov_arg     = params.rt_harmonize_min_anchor_coverage != null ?
                        "--min-anchor-coverage ${params.rt_harmonize_min_anchor_coverage}" : ''
    """
    set -o pipefail

    harmonize_rt.py \\
        --in ${library} \\
        --out ${out_name} \\
        --provenance ${provenance} \\
        --fit-report ${report_name} \\
        --min-anchors ${params.rt_harmonize_min_anchors} \\
        --min-r2 ${params.rt_harmonize_min_r2} \\
        --fallback ${fallback} \\
        ${cov_arg} \\
        ${anchors_arg} \\
        ${pass_flag} \\
        ${full_wnd} \\
        2>&1 | tee harmonize_rt.log

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
      provenance: ${provenance}
      passthrough: ${passthrough}
      emit_full_window: ${emit_full_window}
    EOF
    """

    stub:
    """
    touch ${library.baseName}.irt.${library.extension} \\
          ${library.baseName}.rt_fit.tsv \\
          harmonize_rt.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
