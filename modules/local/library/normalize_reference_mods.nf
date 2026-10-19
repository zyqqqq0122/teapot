process NORMALIZE_REFERENCE_MODS {
    label 'python'
    tag   "${reference_list.baseName}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/prepare_library" },
               mode: params.publish_mode,
               saveAs: { fn -> fn == 'versions.yml' ? 'normalize_reference_mods.versions.yml' : fn }

    input:
    path reference_list
    path library

    output:
    path "${reference_list.baseName}.normalized.${reference_list.extension}", emit: reference_list
    path "normalize_reference_mods.log",                                      emit: log
    path "versions.yml",                                                      emit: versions

    script:
    def out_name = "${reference_list.baseName}.normalized.${reference_list.extension}"
    """
    # bin/normalize_reference_mods.py ${Fingerprint.of("${projectDir}/bin/normalize_reference_mods.py")}: a changed script reruns the step
    set -o pipefail

    normalize_reference_mods.py \\
        --reference ${reference_list} \\
        --library   ${library} \\
        --out       ${out_name} \\
        2>&1 | tee normalize_reference_mods.log

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | sed 's/Python //')
    EOF
    """

    stub:
    def out_name = "${reference_list.baseName}.normalized.${reference_list.extension}"
    """
    cp ${reference_list} ${out_name}
    touch normalize_reference_mods.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
