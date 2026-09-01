process OPENSWATH_ASSAY_GENERATOR {
    label 'openswath'
    tag   "${tsv.baseName}"

    container "${params.openswath_container}"

    publishDir path: { "${params.outdir}/prepare_library_openswath" },
               mode: params.publish_mode

    input:
    path tsv

    output:
    path "${tsv.baseName}.assay.pqp",           emit: pqp
    path "${tsv.baseName}.assay_generator.log", emit: log
    path "versions.yml",                        emit: versions

    script:
    """
    if head -1 ${tsv} | grep -qi "transition_name"; then
        col=\$(head -1 ${tsv} | tr '\t' '\n' | grep -in '^transition_name\$' | cut -d: -f1)
        if [ -n "\$col" ]; then
            dups=\$(awk -F'\t' -v c="\$col" 'NR>1 {print \$c}' ${tsv} | sort | uniq -d)
            if [ -n "\$dups" ]; then
                echo "ASSAY GENERATOR PRE-CHECK FAILED" >&2
                echo "  ${tsv} has \$(printf '%s\n' "\$dups" | wc -l) duplicated transition_name value(s)." >&2
                echo "  OpenSwathAssayGenerator requires them unique and rejects the file." >&2
                echo "  Examples:" >&2
                printf '%s\n' "\$dups" | head -3 | sed 's/^/    /' >&2
                echo "  Make them unique in your library and re-run; the pipeline does not" >&2
                echo "  edit your library for you." >&2
                exit 1
            fi
        fi
    fi

    set -o pipefail

    OpenSwathAssayGenerator \\
        -in  ${tsv} \\
        -out ${tsv.baseName}.assay.pqp \\
        ${Args.flat(params.openswath_assaygen_args)} \\
        2>&1 | tee ${tsv.baseName}.assay_generator.log

    cat <<-EOF > versions.yml
    "${task.process}":
      openms: \$(OpenSwathAssayGenerator --version 2>&1 | head -n1)
    EOF
    """

    stub:
    """
    touch ${tsv.baseName}.assay.pqp ${tsv.baseName}.assay_generator.log
    echo '"${task.process}": {openms: stub}' > versions.yml
    """
}
