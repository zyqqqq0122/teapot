process TRIC_FEATURE_ALIGNMENT {
    label 'tric'
    tag   "tric"

    container 'docker://grosenberger/msproteomicstools:latest'

    publishDir path: { "${params.outdir}/tric" },
               mode: params.publish_mode

    input:
    path tsvs

    output:
    path "feature_alignment.tsv",        emit: long
    path "feature_alignment_matrix.tsv", emit: matrix
    path "tric.log",                     emit: log
    path "versions.yml",                 emit: versions

    script:
    """
    # Fail on a broken pipeline, not just a broken last stage: the command below
    # is piped into `tee`, whose success would otherwise hide the tool's failure.
    set -o pipefail

    feature_alignment.py \\
        --in ${tsvs} \\
        --out feature_alignment.tsv \\
        --out_matrix feature_alignment_matrix.tsv \\
        ${Args.flat(params.tric_args)} \\
        2>&1 | tee tric.log

    cat <<-EOF > versions.yml
    "${task.process}":
      msproteomicstools: \$(feature_alignment.py --help 2>&1 | head -n1 || echo unknown)
    EOF
    """

    stub:
    """
    touch feature_alignment.tsv feature_alignment_matrix.tsv tric.log
    echo '"${task.process}": {msproteomicstools: stub}' > versions.yml
    """
}
