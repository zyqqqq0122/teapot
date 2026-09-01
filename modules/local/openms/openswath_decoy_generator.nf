process OPENSWATH_DECOY_GENERATOR {
    label 'openswath'
    tag   "${pqp.baseName}"

    container "${params.openswath_container}"

    publishDir path: { "${params.outdir}/prepare_library_openswath" },
               mode: params.publish_mode

    input:
    path pqp

    output:
    path "library.pqp",                         emit: pqp
    path "${pqp.baseName}.decoy_generator.log", emit: log
    path "versions.yml",                        emit: versions

    script:
    """
    set -o pipefail

    OpenSwathDecoyGenerator \\
        -in  ${pqp} \\
        -out library.pqp \\
        ${Args.flat(params.openswath_decoygen_args)} \\
        2>&1 | tee ${pqp.baseName}.decoy_generator.log

    if grep -q "Number of target peptides: 0" ${pqp.baseName}.decoy_generator.log; then
        echo "OPENSWATH_DECOY_GENERATOR: the assay library is EMPTY (0 target peptides)." >&2
        echo "  OpenSwathAssayGenerator filtered every precursor out. The usual cause is" >&2
        echo "  -min_transitions in openswath_assaygen_args exceeding the transitions" >&2
        echo "  the library actually carries per precursor." >&2
        exit 1
    fi
    cat <<-EOF > versions.yml
    "${task.process}":
      openms: \$(OpenSwathDecoyGenerator --version 2>&1 | head -n1)
    EOF
    """

    stub:
    """
    touch library.pqp ${pqp.baseName}.decoy_generator.log
    echo '"${task.process}": {openms: stub}' > versions.yml
    """
}
