process OPENSWATH_WORKFLOW {
    label 'openswath'
    tag   "${meta.id}"

    container "${params.openswath_container}"

    publishDir path: { "${params.outdir}/openswath_workflow/${meta.id}" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(mzml), path(reference_list)
    path  library
    path  irt_library, stageAs: 'irt_lib_input/*'
    path  full_window_flag, stageAs: 'full_window/*'
    val   matching_window_only

    output:
    tuple val(meta), path(reference_list), path("${mzml.baseName}.osw"),  emit: osw
    tuple val(meta), path("${meta.id}.openswath.log"),               emit: log
    path  "versions.yml", emit: versions

    script:
    def basename = { f -> f.name.tokenize('/').last() }
    def irt_convert = ''
    def irt_arg     = ''
    if (basename(irt_library) != 'NO_FILE') {
        irt_convert = "TargetedFileConverter -in ${irt_library} -out irt_lib.TraML 2>&1 | tee ${meta.id}.irt_convert.log"
        irt_arg     = '-tr_irt irt_lib.TraML'
    }
    def search_args = Args.flat(params.openswath_search_args)
    def search_args_no_ms1 = (search_args.replaceAll(/-enable_ms1\s+\S+/, '') + ' -enable_ms1 false').trim()
    def force_full_window = ''
    if (basename(full_window_flag) != 'NO_FILE') {
        force_full_window = '-rt_extraction_window -1'
    }
    def warn_full_window = force_full_window ? 'true' : 'false'
    def pinned = (search_args =~ /-matching_window_only\b/)
    def matching_window = (matching_window_only && !pinned) ? '-matching_window_only true' : ''
    """
    set -o pipefail

    if [ "${warn_full_window}" = "true" ]; then
        echo "OPENSWATH_WORKFLOW WARNING: HARMONIZE_RT fell back to full-window extraction for this library. Overriding rt_extraction_window with -1." >&2
    fi

    ${irt_convert}

    openswath() {
        OpenSwathWorkflow \\
            -in     ${mzml} \\
            -tr     ${library} \\
            ${irt_arg} \\
            -out_features ${mzml.baseName}.osw \\
            -threads ${task.cpus} \\
            "\$@" \\
            ${force_full_window} \\
            ${matching_window} \\
            2>&1 | tee ${meta.id}.openswath.log
    }

    if grep -q -m1 -E 'MS:1000511"[^>]*value="1"' ${mzml}; then
        openswath ${search_args}
    else
        echo "OPENSWATH_WORKFLOW: ${mzml} has no MS1 scan; MS1 scoring is off" >&2
        openswath ${search_args_no_ms1}
    fi

    cat <<-EOF > versions.yml
    "${task.process}":
      openms: \$( (OpenSwathWorkflow --helphelp 2>&1 | grep -m1 -i '^Version:') || echo 'Version: unknown')
      tr_irt: ${irt_library.name}
      full_window_override: ${basename(full_window_flag) != 'NO_FILE'}
      matching_window_only: ${matching_window ? 'true' : (pinned ? 'set in openswath_search_args' : 'false')}
    EOF
    """

    stub:
    """
    touch ${mzml.baseName}.osw ${meta.id}.openswath.log
    echo '"${task.process}": {openms: stub}' > versions.yml
    """
}
