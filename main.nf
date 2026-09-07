#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PREPARE_LIBRARY_ENCYCLOPEDIA } from './subworkflows/local/prepare_library_encyclopedia'
include { PREPARE_LIBRARY_OPENSWATH    } from './subworkflows/local/prepare_library_openswath'
include { ENCYCLOPEDIA_DIA             } from './subworkflows/local/encyclopedia_dia'
include { ENCYCLOPEDIA_PRM             } from './subworkflows/local/encyclopedia_prm'
include { OPENSWATH_DIA                } from './subworkflows/local/openswath_dia'
include { OPENSWATH_PRM                } from './subworkflows/local/openswath_prm'
include { MSCONVERT                    } from './modules/local/msconvert/msconvert'


def helpMessage() {
    log.info """
    Usage:
      nextflow run main.nf --input samplesheet.csv --fasta proteins.fasta [options]

    Samplesheet columns (CSV with header):
      sample_id   unique identifier
      file        path to raw MS file: .raw | .mzML | .mzXML | .dia (encyclopedia only)
      mode        DIA | PRM
      reference_list   path to assay .txt/.csv defining the reference peptide
                       set for PRM (used by SPLIT_OPENSWATH_FEATURES / CONTEXT_SEARCH
                       to partition features into reference vs background;
                       ignored for DIA). `reference_list` is accepted as a
                       deprecated alias. May be left empty when a vendor
                       .blib is provided. The pipeline derives it via
                       BLIB_TO_REFERENCE_LIST.
      engine      encyclopedia | openswath  (optional, defaults to encyclopedia)

    Quantification (all engines/modes):
      Base output at results/base_quant/<engine>_<mode>/ has multi-abundance
      columns: abundance_openswath, abundance_encyclopedia, abundance_diathem,
      abundance_tric (when populated), plus abundance_primary chosen by
      --primary_abundance. Non-populated columns are NA.

      --primary_abundance   openswath | encyclopedia | diathem | tric
                            Which column feeds layer-3 calibration.
      --run_diathem         true|false. Off skips diathem; other columns unaffected.
      --diathem_library     Optional .elib/.dlib/.blib for diathem f_prior.
      --min_consistency     LOD gate on diathem-based abs quant (0..1, default 0).

    DIA absolute quantitation (spike-in):
      --heavy_label         List of { residue, unimod } rules. Default SILAC.
      --standard_amounts    TSV/CSV mapping peptide -> known heavy pmol.
    """.stripIndent()
}

workflow {
    if (!params.input || !params.fasta) { helpMessage(); exit 1 }

    def missing_assets = ['NO_FILE', 'unimod_delta_masses.tsv', 'masslist-decoys.jar']
        .findAll { !file("${projectDir}/assets/${it}").exists() }
    if (missing_assets) {
        error "Pre-flight: missing shipped asset(s) in ${projectDir}/assets: " +
              missing_assets.join(', ')
    }

    def renamed = [
        'channel_complete_ion_types':           'complete_channels_ion_types',
        'channel_complete_max_frag_charge':     'complete_channels_max_frag_charge',
        'channel_complete_tolerance_ppm':       'complete_channels_tolerance_ppm',
        'channel_complete_coverage_hard_floor': 'complete_channels_coverage_hard_floor',
        'channel_complete_coverage_warn_below': 'complete_channels_coverage_warn_below',
        'channel_complete_skip_proteins':       'complete_channels_skip_proteins',
    ].findAll { old, _new -> params.containsKey(old) }
    if (renamed) {
        error "Pre-flight: renamed param(s):\n" +
              renamed.collect { o, n -> "  ${o} -> ${n}" }.join('\n')
    }

    def hasPlaceholder
    hasPlaceholder = { x ->
        if (x instanceof CharSequence) return (x =~ /<[\w .\/-]+>/) as Boolean
        if (x instanceof Map)          return x.values().any { hasPlaceholder(it) }
        if (x instanceof Collection)   return x.any { hasPlaceholder(it) }
        return false
    }
    def unfilled = params.findAll { k, v -> hasPlaceholder(v) }
    if (unfilled) {
        error "Pre-flight: unfilled template placeholder(s):\n" +
              unfilled.collect { k, v -> "  ${k} = '${v}'" }.join('\n')
    }

    def rows = file(params.input).splitCsv(header: true)
    def ref_of = { r -> r.reference_list?.trim() ?: '' }
    if (rows && rows.any { it.containsKey('mass_list') } && !rows[0].containsKey('reference_list')) {
        error "samplesheet: the 'mass_list' column was renamed to 'reference_list'."
    }
    def sig  = { r -> [ (r.engine?.trim() ?: 'encyclopedia').toLowerCase(),
                        r.mode.toUpperCase() ] }
    def routes = rows.collect(sig) as Set
    def has_enc     = routes.any { it[0] == 'encyclopedia' }
    def has_osw     = routes.any { it[0] == 'openswath' }
    def has_osw_dia = routes.any { it[0] == 'openswath' && it[1] == 'DIA' }
    def has_osw_prm = routes.any { it[0] == 'openswath' && it[1] == 'PRM' }
    def has_enc_dia = routes.any { it[0] == 'encyclopedia' && it[1] == 'DIA' }
    def has_enc_prm = routes.any { it[0] == 'encyclopedia' && it[1] == 'PRM' }
    def osw_labels = has_osw ? (params.heavy_label ?: []) : []

    def bg_source = params.background_source
    if (params.background_koina && bg_source == 'none') {
        log.warn "background_koina is deprecated; use background_source='koina'"
        bg_source = 'koina'
    }
    def valid_bg_sources = ['none','koina','dlib','blib','elib','gpf']
    if (!(bg_source in valid_bg_sources)) {
        error "background_source must be one of ${valid_bg_sources}, got: '${bg_source}'"
    }

    def enc_prm_uses_derived = has_enc_prm && params.blib && rows.any {
        (it.engine?.trim() ?: 'encyclopedia').toLowerCase() == 'encyclopedia' &&
        it.mode.toUpperCase() == 'PRM' && !ref_of(it)
    }
    def openswath_prm_uses_derived = has_osw_prm && params.openswath_blib && rows.any {
        (it.engine?.trim() ?: 'encyclopedia').toLowerCase() == 'openswath' &&
        it.mode.toUpperCase() == 'PRM' && !ref_of(it)
    }
    def openswath_prm_needs_bg = openswath_prm_uses_derived &&
                           (params.openswath_prm_fdr ?: 'context') != 'pyprophet'
    if ((enc_prm_uses_derived || openswath_prm_needs_bg) && bg_source == 'none') {
        error """
            Pre-flight: PRM sample(s) will use a reference_list derived from the
            vendor .blib (samplesheet reference_list column empty). Without a
            background library the reference covers the entire target library
            and Context's mProphet LDA has zero background targets to train on.

            For the OpenSWATH PRM path you can instead set
              openswath_prm_fdr = 'pyprophet'
            which scores target vs decoy within the search and needs no
            background at all.

            Configure a background source:
              background_source = 'koina'  + background_fasta  (proteome FASTA)
              background_source = 'blib'   + background_library  (background .blib)
              background_source = 'dlib'   + background_library  (background .dlib)
              background_source = 'elib'   + background_library  (pre-built
                                            chromatogram .elib from any
                                            source. Third-party GPF is fine)
              background_source = 'gpf'    + library_sheet + (dlib or
                                            use_koina + background_fasta),
                                            build background from GPF
                                            fractions on the fly. Vendor
                                            blib must NOT be in the seed.
        """.stripIndent().trim()
    }
    def has_vendor_blib = params.blib || params.openswath_blib
    if (has_vendor_blib && params.library_sheet && bg_source != 'gpf') {
        error """
            Pre-flight: a vendor .blib and --library_sheet were both supplied,
            but background_source is '${bg_source}'. GPF fractions would then be
            searched with the VENDOR library as the seed and the resulting
            chromatogram .elib would become the search library, so any heavy
            peptidoform not detected in the GPF runs is silently lost.

            The GPF data is a BACKGROUND:
              --blib / --openswath_blib <vendor.blib>   (target)
              background_source = 'gpf'
              --library_sheet <fractions.csv>
              + a proteome-wide seed: --dlib, or --use_koina + --background_fasta

            (Or drop --library_sheet if you did not mean to build from GPF.)
        """.stripIndent().trim()
    }
    if (has_vendor_blib && params.elib && bg_source != 'elib') {
        error """
            Pre-flight: a vendor .blib and --elib were both supplied, but
            background_source is '${bg_source}'. The .elib would be used as the
            search library and the vendor blib silently discarded, losing its
            heavy channel and its derived reference_list.

            A pre-built GPF chromatogram library is a BACKGROUND. Use:
              --blib / --openswath_blib <vendor.blib>   (target)
              background_source = 'elib'
              --background_library <chromatogram.elib>
        """.stripIndent().trim()
    }

    if (params.openswath_prm_encyclopedia_quant && has_osw_prm &&
        (params.openswath_pqp || params.openswath_tsv)) {
        error """
            Pre-flight: openswath_prm_encyclopedia_quant=true needs an encyclopedia-
            readable library, which the pipeline builds only on the
            --openswath_blib and --use_koina routes. With --openswath_pqp /
            --openswath_tsv there is no .dlib for CONTEXT_SEARCH to search.
            Drop openswath_prm_encyclopedia_quant, or switch to a blib/koina route.
        """.stripIndent().trim()
    }

    if (bg_source == 'koina' && !params.background_fasta) {
        error """
            Pre-flight: background_source='koina' requires --background_fasta.
            Do NOT reuse --fasta if it is a small vendor construct FASTA.
            Predicting a 'background' from vendor constructs yields a
            background composed entirely of vendor targets, which collides
            with the target library and gives zero real background targets.
            Supply a proteome-wide FASTA.
        """.stripIndent().trim()
    }

    def diathemable = { engine, mode -> rows.any { r ->
        (r.engine?.trim() ?: 'encyclopedia').toLowerCase() == engine &&
        r.mode.toUpperCase() == mode &&
        (r.file.toLowerCase() ==~ /.*\.(raw|mzml|mzxml)$/)
    } }
    def diathem_ok_enc_dia = params.run_diathem && diathemable('encyclopedia', 'DIA')
    def diathem_ok_enc_prm = params.run_diathem && diathemable('encyclopedia', 'PRM')
    def diathem_ok_osw_dia = params.run_diathem && has_osw_dia
    def diathem_ok_osw_prm = params.run_diathem && has_osw_prm

    def runs_by_route = rows
        .groupBy { r -> [ (r.engine?.trim() ?: 'encyclopedia').toLowerCase(),
                          r.mode.toUpperCase() ] }
        .collectEntries { k, v -> [ (k): v.size() ] }
    runs_by_route.each { route, n ->
        if (n < params.quant_min_runs) {
            def msg = "quantification cohort for ${route[0]}/${route[1]} is ${n} " +
                      "run(s), below quant_min_runs=${params.quant_min_runs}. " +
                      "Cross-run quantification (EncyclopeDIA RT alignment, " +
                      "diathem consistency, pyprophet, TRIC) is unreliable at " +
                      "this cohort size; identifications are unaffected."
            if (params.quant_require_min_runs) {
                error """
                    Pre-flight: ${msg}

                    Either supply more runs, lower quant_min_runs, or set
                    quant_require_min_runs = false to downgrade this to a warning.
                """.stripIndent().trim()
            }
            log.warn "Pre-flight: ${msg} Treat abundances as within-run values."
        }
    }

    def sample_map_str = rows.collect { r ->
        def f = file(r.file)
        def stem = f.getName().replaceFirst(/\.(raw|mzML|mzml|mzXML|mzxml|dia)$/, '')
        "${stem}\t${r.sample_id}"
    }.join('\n') + '\n'
    def sample_map_file = file("${workDir}/sample_map.tsv")
    sample_map_file.parent.mkdirs()
    sample_map_file.text = sample_map_str
    def sample_map_ch = Channel.value(sample_map_file)

    samples = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def engine = (row.engine?.trim() ?: 'encyclopedia').toLowerCase()
            if (!(engine in ['encyclopedia', 'openswath'])) {
                error "samplesheet: unknown engine '${engine}' for ${row.sample_id}"
            }
            def meta = [id: row.sample_id, mode: row.mode.toUpperCase(), engine: engine]
            def f    = file(row.file)
            def rl_val = ref_of(row)
            def rl   = rl_val ? file(rl_val) : file("${projectDir}/assets/NO_FILE")
            if (engine == 'openswath' && f.getName().toLowerCase().endsWith('.dia')) {
                error ".dia is encyclopedia-specific; sample '${meta.id}' " +
                      "uses engine=openswath. Provide the .raw or a .mzML instead."
            }
            tuple(meta, f, rl)
        }

    fasta            = file(params.fasta)
    no_file          = file("${projectDir}/assets/NO_FILE")
    standard_amts_ch = params.standard_amounts ? file(params.standard_amounts) : no_file
    diathem_lib_ch   = params.diathem_library ? file(params.diathem_library) : no_file

    def VENDOR_EXT = ['.raw', '.wiff', '.wiff2', '.d', '.lcd']
    def is_vendor  = { f -> VENDOR_EXT.any { f.getName().toLowerCase().endsWith(it) } }

    def sidecars_of = { f ->
        ['.scan', '2.scan', '.timeseries.data'].collect { file("${f}${it}") }
                                               .findAll { it.exists() }
    }

    ext_split = samples.branch { meta, f, _ml ->
        raw:   is_vendor(f)
        ready: true
    }

    MSCONVERT(ext_split.raw.map { meta, f, rl -> tuple(meta, f, sidecars_of(f), rl) })

    converted = MSCONVERT.out.mzml.flatMap { meta, mzml, rl ->
        def files = (mzml instanceof List) ? mzml : [mzml]
        if (files.size() == 1) return [ tuple(meta, files[0], rl) ]
        log.warn "MSCONVERT: ${meta.id} is a multi-sample acquisition -> " +
                 "${files.size()} runs: ${files.collect{ it.baseName }.join(', ')}"
        files.collect { f -> tuple(meta + [id: f.baseName], f, rl) }
    }

    normalized = ext_split.ready.mix(converted)

    by_route = normalized.branch { meta, _f, _ml ->
        enc_dia: meta.engine == 'encyclopedia' && meta.mode == 'DIA'
        enc_prm: meta.engine == 'encyclopedia' && meta.mode == 'PRM'
        osw_dia: meta.engine == 'openswath'    && meta.mode == 'DIA'
        osw_prm: meta.engine == 'openswath'    && meta.mode == 'PRM'
    }

    def enc_prm_row = has_enc_prm
        ? rows.find {
            (it.engine?.trim() ?: 'encyclopedia').toLowerCase() == 'encyclopedia' &&
            it.mode.toUpperCase() == 'PRM' }
        : null
    def _enc_prm_ref = enc_prm_row ? ref_of(enc_prm_row) : ''
    def enc_prm_sample_rl = _enc_prm_ref ? file(_enc_prm_ref) : null

    def enc_dia_row = has_enc_dia
        ? rows.find {
            (it.engine?.trim() ?: 'encyclopedia').toLowerCase() == 'encyclopedia' &&
            it.mode.toUpperCase() == 'DIA' && ref_of(it) }
        : null
    def _enc_dia_ref = enc_dia_row ? ref_of(enc_dia_row) : ''
    def enc_dia_sample_rl = _enc_dia_ref ? file(_enc_dia_ref) : null

    def enc_labels = has_enc ? (params.heavy_label ?: []) : []

    if (has_enc) {
        PREPARE_LIBRARY_ENCYCLOPEDIA(fasta, enc_labels, bg_source,
                        params.background_library, params.background_fasta,
                        params.background_min_targets, no_file)

        def enc_prm_samples
        def enc_prm_targets
        if (has_enc_prm) {
            if (enc_prm_sample_rl != null) {
                enc_prm_samples = by_route.enc_prm
                enc_prm_targets = Channel.value(enc_prm_sample_rl)
            } else if (params.blib) {
                enc_prm_samples = by_route.enc_prm
                    .combine(PREPARE_LIBRARY_ENCYCLOPEDIA.out.reference_list_derived)
                    .map { meta, f, _orig_rl, derived -> tuple(meta, f, derived) }
                enc_prm_targets = PREPARE_LIBRARY_ENCYCLOPEDIA.out.reference_list_derived
            } else {
                error "encyclopedia+PRM samples need a reference_list, either " +
                      "in the samplesheet's reference_list column, or via " +
                      "--blib for automatic " +
                      "derivation from a vendor .blib."
            }
        } else {
            enc_prm_samples = by_route.enc_prm
            enc_prm_targets = no_file
        }

        def enc_dia_targets = no_file
        if (has_enc_dia) {
            if (enc_dia_sample_rl != null)      enc_dia_targets = Channel.value(enc_dia_sample_rl)
            else if (params.blib)               enc_dia_targets = PREPARE_LIBRARY_ENCYCLOPEDIA.out.reference_list_derived
        }

        ENCYCLOPEDIA_DIA(by_route.enc_dia, PREPARE_LIBRARY_ENCYCLOPEDIA.out.library, fasta,
                     diathem_lib_ch, standard_amts_ch, params.heavy_label,
                     sample_map_ch, enc_dia_targets, no_file, diathem_ok_enc_dia)
        ENCYCLOPEDIA_PRM(enc_prm_samples, PREPARE_LIBRARY_ENCYCLOPEDIA.out.library, fasta,
                     diathem_lib_ch, standard_amts_ch, params.heavy_label,
                     sample_map_ch, enc_prm_targets, no_file, diathem_ok_enc_prm)
    }

    if (has_osw) {
        def openswath_prm_refs = has_osw_prm ? rows.findAll {
            (it.engine?.trim() ?: 'encyclopedia').toLowerCase() == 'openswath' &&
            it.mode.toUpperCase() == 'PRM' && ref_of(it)
        }.collect { ref_of(it) }.unique() : []
        if (openswath_prm_refs.size() > 1) {
            error """
                main.nf: OSW-PRM samples in the samplesheet reference more than
                one distinct reference_list file:
                  ${openswath_prm_refs.join('\n                  ')}
                The pipeline derives iRT anchors from the reference_list and
                applies them library-wide. Using one sample's list to calibrate
                another sample's data would silently cross-calibrate.
                Split into separate pipeline runs, one per reference_list,
                or make every OSW-PRM row reference the same file.
                """.stripIndent().trim()
        }
        def openswath_reflist_for_anchors = openswath_prm_refs ?
                                          file(openswath_prm_refs[0]) : no_file

        PREPARE_LIBRARY_OPENSWATH(fasta, osw_labels, bg_source,
                            params.background_library, params.background_fasta,
                            params.background_min_targets,
                            openswath_reflist_for_anchors, no_file)

        if (has_osw_prm) {
            def openswath_prm_row = rows.find {
                (it.engine?.trim() ?: 'encyclopedia').toLowerCase() == 'openswath' &&
                it.mode.toUpperCase() == 'PRM'
            }
            def _openswath_ref = ref_of(openswath_prm_row)
            def sample_rl = _openswath_ref ? file(_openswath_ref) : null

            def openswath_prm_samples
            def openswath_prm_targets
            if (sample_rl != null) {
                openswath_prm_samples = by_route.osw_prm
                openswath_prm_targets = Channel.value(sample_rl)
            } else if (params.openswath_blib) {
                openswath_prm_samples = by_route.osw_prm
                    .combine(PREPARE_LIBRARY_OPENSWATH.out.reference_list_derived)
                    .map { meta, f, _orig_rl, derived -> tuple(meta, f, derived) }
                openswath_prm_targets = PREPARE_LIBRARY_OPENSWATH.out.reference_list_derived
            } else {
                error "openswath+PRM samples need a reference_list, either " +
                      "in the samplesheet's reference_list column, or via " +
                      "--openswath_blib for automatic " +
                      "derivation from a vendor .blib."
            }

            OPENSWATH_PRM(openswath_prm_samples,
                                   PREPARE_LIBRARY_OPENSWATH.out.library,
                                   PREPARE_LIBRARY_OPENSWATH.out.irt_library,
                                   PREPARE_LIBRARY_OPENSWATH.out.full_window_flag,
                                   PREPARE_LIBRARY_OPENSWATH.out.encyclopedia_library,
                                   fasta,
                                   diathem_lib_ch, standard_amts_ch,
                                   params.heavy_label, sample_map_ch,
                                   openswath_prm_targets, no_file, diathem_ok_osw_prm)
        }
        if (has_osw_dia) {
            OPENSWATH_DIA(by_route.osw_dia,
                                   PREPARE_LIBRARY_OPENSWATH.out.library,
                                   PREPARE_LIBRARY_OPENSWATH.out.irt_library,
                                   PREPARE_LIBRARY_OPENSWATH.out.full_window_flag,
                                   fasta,
                                   diathem_lib_ch, standard_amts_ch,
                                   params.heavy_label, sample_map_ch,
                                   no_file, diathem_ok_osw_dia)
        }
    }
}
