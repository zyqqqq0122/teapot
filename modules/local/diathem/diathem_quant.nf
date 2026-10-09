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
    val   rt_window

    output:
    tuple val(mode), path("${mode}.diathem_quant.tsv"),  emit: quant
    tuple val(mode), path("${mode}.sample_map.tsv"),     emit: sample_map
    path  "${mode}.diathem.log",                          emit: log
    path  "versions.yml",                                 emit: versions

    script:
    def def_tol = mode == 'PRM' ? '0.5' : '12.5'
    def tol_cmd = params.diathem_precursor_tol != null
        ? "tol=${params.diathem_precursor_tol}"
        : "tol=\$(diathem_quant.py --targets ${targets} --fallback ${def_tol} mzml_dir/*.mzML 2>> ${mode}.diathem.log)"
    def rtw  = params.diathem_rt_window     != null ? "${params.diathem_rt_window}"     : "${rt_window}"
    def fth   = params.diathem_fragment_tol_th != null ? "--fragment-tol-th ${params.diathem_fragment_tol_th}" : ''
    def fppm  = (params.diathem_fragment_ppm != null && params.diathem_fragment_tol_th == null) ? "--fragment-ppm ${params.diathem_fragment_ppm}" : ''
    def lppm  = params.diathem_library_ppm  != null ? "--library-ppm ${params.diathem_library_ppm}"   : ''
    def extra = Args.flat(params.diathem_args)
    def want_weight = !(extra =~ /--weight\b/)
    def want_top = (params.diathem_max_fragments as int) > 0 && !(extra =~ /--max-fragments\b/)
    def want_bl  = params.diathem_baseline && !(extra =~ /--background\b/)
    def flags = "--rt-window ${rtw} ${fth} ${fppm} ${lppm}".replaceAll(/\s+/, ' ').trim()
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
    if [ "\$(ls mzml_dir | wc -l)" -ne "\$(ls mzml_in | wc -l)" ]; then
        echo "DIATHEM_QUANT: \$(ls mzml_in | wc -l) mzML staged, \$(ls mzml_dir | wc -l) matched to a sample in ${sample_map}" >&2
        exit 1
    fi

    : > ${mode}.diathem.log
    ${tol_cmd}
    dq_help=\$(diathem quant --help 2>/dev/null || true)
    bench=''
    if [ "${want_weight}" = "true" ]; then
        if echo "\$dq_help" | grep -o -m1 -- '--weight {[^}]*}' | tr '{,}' '   ' | grep -qw -- '${params.diathem_weight}'; then
            bench="\$bench --weight ${params.diathem_weight}"
        else
            echo "diathem_weight: this diathem has no ${params.diathem_weight} weighting; using uniform" >> ${mode}.diathem.log
            bench="\$bench --weight uniform"
        fi
    fi
    if [ "${want_top}" = "true" ]; then
        if echo "\$dq_help" | grep -q -- '--max-fragments'; then
            bench="\$bench --max-fragments ${params.diathem_max_fragments} --min-fragment-mz 200"
        else
            echo "diathem_max_fragments: this diathem has no --max-fragments; fitting every fragment" >> ${mode}.diathem.log
        fi
    fi
    if [ "${want_bl}" = "true" ]; then
        if echo "\$dq_help" | grep -q -- '--background'; then
            bench="\$bench --background"
        else
            echo "diathem_baseline: this diathem has no --background; fitting without a baseline" >> ${mode}.diathem.log
        fi
    fi
    cat ${mode}.diathem.log >&2

    diathem quant \\
        --targets  ${targets} \\
        --mzml-dir mzml_dir/ \\
        ${lib_arg} \\
        --output   ${mode}.diathem_quant.tsv \\
        --precursor-tol \$tol ${flags} \$bench \\
        ${extra} \\
        2>&1 | tee -a ${mode}.diathem.log

    cp ${sample_map} ${mode}.sample_map.tsv

    cat <<-EOF > versions.yml
    "${task.process}":
      diathem: \$(python3 -c 'import importlib.metadata as m; print(m.version("diathem"))' 2> /dev/null || echo unknown)
      container: ${params.diathem_container}
      mode:    ${mode}
      flags:   '--precursor-tol \$tol ${flags} \$bench ${extra}'
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
