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
    path  nested_flag, stageAs: 'nested_windows/*'

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
    def force_full_window = ''
    if (basename(full_window_flag) != 'NO_FILE') {
        force_full_window = '-rt_extraction_window -1'
    }
    def warn_full_window = force_full_window ? 'true' : 'false'
    def matching_window = (basename(nested_flag) != 'NO_FILE') ? '-matching_window_only true' : ''
    def warn_matching_window = matching_window ? 'true' : 'false'
    """
    # Fail on a broken pipeline, not just a broken last stage: every command
    # below is piped into `tee`, whose success would otherwise hide theirs.
    set -o pipefail

    if [ "${warn_full_window}" = "true" ]; then
        echo "OPENSWATH_WORKFLOW WARNING: HARMONIZE_RT fell back to full-window extraction for this library. Overriding rt_extraction_window with -1." >&2
    fi

    ${irt_convert}

    OpenSwathWorkflow \\
        -in     ${mzml} \\
        -tr     ${library} \\
        ${irt_arg} \\
        -out_features ${mzml.baseName}.osw \\
        -threads ${task.cpus} \\
        ${search_args} \\
        ${force_full_window} \\
        2>&1 | tee ${meta.id}.openswath.log

    cat <<-EOF > versions.yml
    "${task.process}":
      openms: \$( (OpenSwathWorkflow --helphelp 2>&1 | grep -m1 -i '^Version:') || echo 'Version: unknown')
      tr_irt: ${irt_library.name}
      full_window_override: ${basename(full_window_flag) != 'NO_FILE'}
    EOF
    """

    stub:
    """
    touch ${mzml.baseName}.osw ${meta.id}.openswath.log
    echo '"${task.process}": {openms: stub}' > versions.yml
    """
}
