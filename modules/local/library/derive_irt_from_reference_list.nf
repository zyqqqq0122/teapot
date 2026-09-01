process DERIVE_IRT_FROM_REFERENCE_LIST {
    label 'python'
    tag   "${reference_list.baseName}"

    container "${params.python_container}"

    publishDir path: { "${params.outdir}/derive_irt_from_reference_list" },
               mode: params.publish_mode

    input:
    path reference_list
    path koina_dlib

    output:
    path "${reference_list.baseName}.reflist_anchors.tsv",  emit: anchors
    path "${reference_list.baseName}.reflist_irt.tsv",      emit: osw_tsv
    path "derive_irt_from_reference_list.log",              emit: log
    path "versions.yml",                                    emit: versions

    script:
    """
    set -o pipefail

    derive_irt_from_reference_list.py \\
        --reference-list ${reference_list} \\
        --koina-dlib     ${koina_dlib} \\
        --anchors-tsv    ${reference_list.baseName}.reflist_anchors.tsv \\
        --osw-tsv        ${reference_list.baseName}.reflist_irt.tsv \\
        --top-transitions ${params.irt_top_transitions} \\
        2>&1 | tee derive_irt_from_reference_list.log

    cat <<-EOF > versions.yml
    "${task.process}":
      python: \$(python3 --version 2>&1 | awk '{print \$2}')
    EOF
    """

    stub:
    """
    touch ${reference_list.baseName}.reflist_anchors.tsv \\
          ${reference_list.baseName}.reflist_irt.tsv \\
          derive_irt_from_reference_list.log
    echo '"${task.process}": {python: stub}' > versions.yml
    """
}
