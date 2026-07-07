process MSCONVERT {
    label 'msconvert'
    tag   "${meta.id}"

    container        'docker://proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses'
    containerOptions '-B `mktemp -d /dev/shm/wineXXX`:/mywineprefix --writable-tmpfs'

    publishDir path: { "${params.outdir}/msconvert/${meta.id}" },
               mode: params.publish_mode

    input:
    tuple val(meta), path(raw), path(mass_list)

    output:
    tuple val(meta), path("${raw.baseName}.mzML"), path(mass_list), emit: mzml
    tuple val(meta), path("${meta.id}.msconvert.log"),              emit: log
    path  "versions.yml",                                            emit: versions

    script:
    """
    set -euo pipefail
    ulimit -s 8192

    mywine msconvert \\
        ${raw} \\
        --mzML \\
        --verbose \\
        ${params.peakpicking} ${params.otherfilters} \\
        2>&1 | tee ${meta.id}.msconvert.log

    cat <<-EOF > versions.yml
    "${task.process}":
      container: proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses
    EOF
    """

    stub:
    """
    touch ${raw.baseName}.mzML ${meta.id}.msconvert.log
    echo '"${task.process}": {msconvert: stub}' > versions.yml
    """
}
