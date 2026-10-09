process MSCONVERT {
    label 'msconvert'
    tag   "${meta.id}"

    container        "${params.msconvert_container}"
    containerOptions "-B `mktemp -d /dev/shm/wineXXX`:/mywineprefix -B ${projectDir}/assets/wine_cpus_online:/sys/devices/system/cpu/online:ro"

    publishDir path: { "${params.outdir}/msconvert/${meta.id}" },
               mode: params.publish_mode,
               saveAs: { fn -> fn ==~ /.*\.(mzML|log)$/ ? fn : null }

    input:
    tuple val(meta), path(raw), path(sidecars), path(reference_list)

    output:
    tuple val(meta), path("*.mzML"), path(reference_list), emit: mzml
    tuple val(meta), path("${meta.id}.msconvert.log"),              emit: log
    path  "versions.yml",                                            emit: versions

    script:
    """
    set -euo pipefail
    ulimit -s 8192
    trap 'WINEPREFIX=/mywineprefix/ wineserver -k 2> /dev/null || true; find /mywineprefix -mindepth 1 -delete 2> /dev/null || true' EXIT

    mywine msconvert \\
        ${raw} \\
        --mzML \\
        --verbose \\
        ${params.peakpicking} ${params.otherfilters} \\
        2>&1 | tee ${meta.id}.msconvert.log

    cat <<-EOF > versions.yml
    "${task.process}":
      container: ${params.msconvert_container}
    EOF
    """

    stub:
    """
    touch ${raw.baseName}.mzML ${meta.id}.msconvert.log
    echo '"${task.process}": {msconvert: stub}' > versions.yml
    """
}
