process SEARCH_DIA {
    label 'encyclopedia'
    tag   "${meta.id}"

    publishDir path: { "${params.outdir}/search_dia/${meta.id}" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(dia), path(library), path(windows)
    path  fasta

    output:
    tuple val(meta), path("${dia.name}.features.txt"),                              emit: features
    tuple val(meta), path("${dia.name}.encyclopedia.txt"),       optional: true,     emit: targets
    tuple val(meta), path("${dia.name}.encyclopedia.decoy.txt"), optional: true,     emit: decoys
    tuple val(meta), path("${dia.name}.elib"),                    optional: true,     emit: elib
    tuple val(meta), path("${meta.id}.search.log"),                                       emit: log
    path  "versions.yml", emit: versions

    script:
    def jar_fp = Fingerprint.of(params.encyclopedia_jar)
    """
    set -o pipefail

    echo "==> encyclopedia jar: ${params.encyclopedia_jar} [${jar_fp}]"

    java -Xmx${params.java_mem} -cp ${params.encyclopedia_jar} \\
        edu.washington.gs.maccoss.encyclopedia.Encyclopedia \\
        -i ${dia} \\
        -l ${library} \\
        -f ${fasta} \\
        ${Args.flat(params.dia_search_args)} \\
        ${windows.name != 'NO_FILE' ? "-apexWindows ${windows}" : ''} \\
        2>&1 | tee ${meta.id}.search.log

    if [ "${windows.name}" != "NO_FILE" ] && ! grep -q "Apex windows: read" ${meta.id}.search.log; then
        echo "ERROR: ${meta.id} has retention-time windows (rt_alignment), but the EncyclopeDIA in use does" >&2
        echo "not read them (-apexWindows needs contextbook 6.10.19-teapot or later). Update" >&2
        echo "encyclopedia_container, or set rt_alignment: off." >&2
        exit 1
    fi

    cat <<-EOF > versions.yml
    "${task.process}":
      encyclopedia: \$(java -jar ${params.encyclopedia_jar} -version 2>&1 | head -n1)
    EOF
    """

    stub:
    """
    touch ${dia.name}.features.txt \\
          ${dia.name}.encyclopedia.txt \\
          ${dia.name}.encyclopedia.decoy.txt \\
          ${dia.name}.elib \\
          ${meta.id}.search.log
    echo '"${task.process}": {encyclopedia: stub}' > versions.yml
    """
}
