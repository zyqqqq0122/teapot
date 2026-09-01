process DIATHEM_QUANT {
    label 'diathem'
    tag   "${mode}"

    container "${params.diathem_container}"

    publishDir path: { "${params.outdir}/diathem_quant/${mode}" },
               mode: params.publish_mode

    input:
    val   mode
    path  mzmls, stageAs: 'mzml_in/*'
    path  targets
    path  library, stageAs: 'prior_lib/*'
    path  sample_map

    output:
    tuple val(mode), path("${mode}.diathem_quant.tsv"),  emit: quant
    tuple val(mode), path("${mode}.sample_map.tsv"),     emit: sample_map
    path  "${mode}.diathem.log",                          emit: log
    path  "versions.yml",                                 emit: versions

    script:
    def def_tol = mode == 'PRM' ? '0.5' : '12.5'
    def def_rtw = mode == 'PRM' ? '60'  : '0'
    def tol  = params.diathem_precursor_tol != null ? "${params.diathem_precursor_tol}" : def_tol
    def rtw  = params.diathem_rt_window     != null ? "${params.diathem_rt_window}"     : def_rtw
    def fth   = params.diathem_fragment_tol_th != null ? "--fragment-tol-th ${params.diathem_fragment_tol_th}" : ''
    def fppm  = (params.diathem_fragment_ppm != null && params.diathem_fragment_tol_th == null) ? "--fragment-ppm ${params.diathem_fragment_ppm}" : ''
    def lppm  = params.diathem_library_ppm  != null ? "--library-ppm ${params.diathem_library_ppm}"   : ''
    def flags = "--precursor-tol ${tol} --rt-window ${rtw} ${fth} ${fppm} ${lppm}".replaceAll(/\s+/, ' ').trim()
    def srcenv = params.diathem_src ? "export PYTHONPATH=${params.diathem_src}:\${PYTHONPATH:-}" : ':'
    def lib_arg = library.name.tokenize('/').last() == 'NO_FILE' ? '' : "--library ${library}"
    """
    set -o pipefail

    ${srcenv}
    mkdir -p mzml_dir
    while IFS=\$'\\t' read -r stem sid; do
        for ext in mzML mzml; do
            src=mzml_in/\${stem}.\${ext}
            [ -f "\$src" ] && ln -sf "\$(readlink -f "\$src")" mzml_dir/\${sid}.mzML
        done
    done < ${sample_map}

    diathem quant \\
        --targets  ${targets} \\
        --mzml-dir mzml_dir/ \\
        ${lib_arg} \\
        --output   ${mode}.diathem_quant.tsv \\
        ${flags} \\
        ${Args.flat(params.diathem_args)} \\
        2>&1 | tee ${mode}.diathem.log

    cp ${sample_map} ${mode}.sample_map.tsv

    cat <<-EOF > versions.yml
    "${task.process}":
      diathem: \$(diathem --version 2>&1 | tail -n1 || echo unknown)
      mode:    ${mode}
      flags:   '${flags}'
      library: '${library.name}'
      src:     '${params.diathem_src ?: '(image)'}'
    EOF
    """

    stub:
    """
    touch ${mode}.diathem_quant.tsv ${mode}.diathem.log
    cp ${sample_map} ${mode}.sample_map.tsv
    echo '"${task.process}": {diathem: stub, mode: ${mode}}' > versions.yml
    """
}
