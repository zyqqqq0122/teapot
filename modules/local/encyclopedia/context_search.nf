process CONTEXT_SEARCH {
    label 'encyclopedia'
    tag   "${meta.id}"

    publishDir path: { "${params.outdir}/context_search/${meta.id}" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(prm), path(reference_list)
    path  library
    path  fasta

    output:
    tuple val(meta), path("${prm.baseName}_reference.features.txt"),               emit: ref_features
    tuple val(meta), path("${prm.baseName}_background.features.txt"),              emit: bg_features
    tuple val(meta), path("${prm.baseName}_reference.features.pep.output.txt"),    emit: ref_targets
    tuple val(meta), path("${prm.baseName}_reference.features.pep.decoy.txt"),     emit: ref_decoys
    tuple val(meta), path("${prm.baseName}_background.features.pep.output.txt"),   emit: bg_targets
    tuple val(meta), path("${prm.baseName}_background.features.pep.decoy.txt"),    emit: bg_decoys
    tuple val(meta), path("${prm.name}*", includeInputs: true),                    emit: search_products
    tuple val(meta), path("plots"), optional: true,                                emit: plots
    tuple val(meta), path("${meta.id}.context.log"),                                emit: log
    path  "versions.yml", emit: versions

    script:
    def jar_fp = Fingerprint.of(params.encyclopedia_jar)
    """
    set -o pipefail

    echo "==> encyclopedia jar: ${params.encyclopedia_jar} [${jar_fp}]"

    if [ ! -s "${reference_list}" ] || [ "${reference_list.name}" = "NO_FILE" ]; then
        echo "ERROR: PRM sample ${meta.id} requires a reference_list in the samplesheet." >&2
        exit 1
    fi

    java -Xmx${params.java_mem} -cp ${params.encyclopedia_jar} \\
        edu.washington.gs.maccoss.encyclopedia.context.ContextMProphetExecutor \\
        -i        ${prm} \\
        -l        ${library} \\
        -f        ${fasta} \\
        -massList ${reference_list} \\
        -fdr      ${params.fdr} \\
        -seed     1 \\
        -plotsdir plots \\
        ${Args.flat(params.context_search_args)} \\
        2>&1 | tee ${meta.id}.context.log

    cp ${prm.baseName}_reference.features.pep.output.txt ${prm.name}.encyclopedia.txt
    cp ${prm.baseName}_reference.features.pep.decoy.txt  ${prm.name}.encyclopedia.decoy.txt

    cat <<-EOF > versions.yml
    "${task.process}":
      encyclopedia: \$(java -jar ${params.encyclopedia_jar} -version 2>&1 | head -n1)
      fdr: ${params.fdr}
    EOF
    """

    stub:
    """
    touch ${prm.baseName}_reference.features.txt \\
          ${prm.baseName}_background.features.txt \\
          ${prm.baseName}_reference.features.pep.output.txt \\
          ${prm.baseName}_reference.features.pep.decoy.txt \\
          ${prm.baseName}_background.features.pep.output.txt \\
          ${prm.baseName}_background.features.pep.decoy.txt \\
          ${meta.id}.context.log \\
          ${prm.name}.features.txt \\
          ${prm.name}.encyclopedia.txt \\
          ${prm.name}.encyclopedia.decoy.txt
    mkdir -p plots
    echo '"${task.process}": {encyclopedia: stub}' > versions.yml
    """
}
