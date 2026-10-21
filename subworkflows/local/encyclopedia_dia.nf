include { SEARCH_DIA                   } from '../../modules/local/encyclopedia/search_dia'
include { ALIGN_RETENTION_TIMES        } from '../../modules/local/runs/align_retention_times'
include { ENCYCLOPEDIA_LIBEXPORT       } from '../../modules/local/encyclopedia/encyclopedia_libexport'
include { ASSERT_QUANT_PRODUCTIVE      } from '../../modules/local/quant/assert_quant_productive'
include { DIATHEM_QUANT                } from '../../modules/local/diathem/diathem_quant'
include { SELECT_DIATHEM_TARGETS       } from '../../modules/local/diathem/select_diathem_targets'
include { STRIP_REFERENCE_DECOYS      } from '../../modules/local/library/strip_reference_decoys'
include { MERGE_QUANT_ENCYCLOPEDIA_DIA } from '../../modules/local/quant/merge_quant_encyclopedia_dia'
include { DERIVE_PAIR_RT_TOLERANCE     } from '../../modules/local/quant/derive_pair_rt_tolerance'
include { DERIVE_PAIR_MIN_COELUTION    } from '../../modules/local/quant/derive_pair_min_coelution'
include { DERIVE_SEARCH_TOLERANCE      } from '../../modules/local/runs/derive_search_tolerance'
include { FINALIZE_QUANT               } from '../../modules/local/quant/finalize_quant'
include { QUANTIFY_HEAVY_LIGHT            } from '../../modules/local/quant/quantify_heavy_light'
include { QUANTIFY_STANDARD_PEAK       } from '../../modules/local/quant/quantify_standard_peak'

workflow ENCYCLOPEDIA_DIA {
    take:
    samples
    library
    fasta
    diathem_library
    standard_amounts
    heavy_label
    reference_list
    has_reference_list
    no_file
    diathem_ok
    located
    acquired
    fragment_errors
    rt_alignment
    enc_quant

    main:
    sample_map = samples
        .map { meta, f, _ml -> "${f.baseName}\t${meta.id}" }
        .collectFile(name: 'sample_map.enc_dia.tsv', newLine: true, sort: true)
        .first()

    def sil = (heavy_label as List)?.size() > 0
    def warn_ch = Channel.empty()

    def runs_in
    if (rt_alignment && sil && has_reference_list) {
        def ids = samples.map { meta, _f, _ml -> meta.id }
            .collectFile(name: 'samples.enc_dia.txt', newLine: true, sort: true)
        def mine = samples.map { meta, _f, _ml -> [meta.id] }
            .join(located.map { meta, t -> [meta.id, t] })
            .map { _id, t -> t }
        ALIGN_RETENTION_TIMES(ids, mine.collect().ifEmpty([]), acquired.collect().ifEmpty([]),
                              library, reference_list, heavy_label, 'enc_dia')
        warn_ch = warn_ch.mix(ALIGN_RETENTION_TIMES.out.warnings.flatMap { Warnings.inFile(it, 'enc_dia') })
        def windows = ALIGN_RETENTION_TIMES.out.windows.flatten()
            .map { w -> [w.name - '.apex_windows.tsv', w] }
        runs_in = samples.map { meta, f, _ml -> [meta.id, meta, f] }
            .join(windows, remainder: true)
            .map { _id, meta, f, w -> tuple(meta, f, w ?: no_file) }
    } else {
        runs_in = samples.map { meta, f, _ml -> tuple(meta, f, no_file) }
    }
    DERIVE_SEARCH_TOLERANCE(fragment_errors.map { meta, f -> tuple(meta.id, f) }
                                .join(samples.map { meta, _f, _ml -> tuple(meta.id, true) })
                                .map { _id, f, _s -> f }.collect().ifEmpty([no_file]),
                            Args.fragmentTolerance(params.dia_search_args), 'enc_dia')
    def tol_ch = DERIVE_SEARCH_TOLERANCE.out.tolerance.combine(DERIVE_SEARCH_TOLERANCE.out.runs)
        .map { t, r -> Args.readTolerance(t, r) }
        .first()
    def args_of = { Map t, String sid -> t.runs.contains(sid) ? Args.withBatchTolerance(params.dia_search_args, t) : Args.flat(params.dia_search_args) }
    def batch_args = tol_ch.map { t -> Args.withBatchTolerance(params.dia_search_args, t) }
    tol_ch.subscribe { t -> log.info Warnings.fragmentTolerance(t, 'enc_dia', params.dia_search_args) }
    SEARCH_DIA(runs_in.combine(library).combine(tol_ch).map { meta, dia, w, lib, t -> tuple(meta, dia, lib, w, args_of(t, meta.id)) }, fasta)

    def elib_ch     = Channel.value(no_file)
    def peptides_ch = Channel.value(no_file)
    if (enc_quant) {
        pooled = SEARCH_DIA.out.features.map { _m, f -> f }
            .mix(SEARCH_DIA.out.targets.map { _m, t -> t })
            .mix(SEARCH_DIA.out.decoys .map { _m, d -> d })
            .mix(SEARCH_DIA.out.elib   .map { _m, e -> e })
            .mix(samples.map { _m, f, _ml -> f })
            .collect(sort: true)
        ENCYCLOPEDIA_LIBEXPORT(pooled, 'analysis', library, fasta, batch_args,
                               has_reference_list)

        ASSERT_QUANT_PRODUCTIVE(ENCYCLOPEDIA_LIBEXPORT.out.elib, 'encyclopedia_dia')
        warn_ch = warn_ch.mix(ASSERT_QUANT_PRODUCTIVE.out.warnings.flatMap { Warnings.inFile(it, 'enc_dia') })
        elib_ch = ENCYCLOPEDIA_LIBEXPORT.out.elib
        peptides_ch = ENCYCLOPEDIA_LIBEXPORT.out.peptides
            .ifEmpty {
                log.warn "ENCYCLOPEDIA_DIA: libexport wrote no analysis.elib.peptides.txt " +
                         "(single search, no cross-run alignment). Abundance comes " +
                         "from the .elib peptidequants table, which is the primary " +
                         "source in either case."
                return file("${projectDir}/assets/NO_FILE")
            }
    }

    enc_targets_ch = SEARCH_DIA.out.targets.map { _m, t -> t }.collect(sort: true)

    def diathem_tsv_ch
    if (diathem_ok) {
        mzmls_ch = samples
            .filter { _m, f, _ml -> f.getName().toLowerCase().endsWith('.mzml') }
            .map    { _m, f, _ml -> f }
            .collect(sort: true)
        def dia_prior = (diathem_library.name != 'NO_FILE') ? diathem_library : library
        def dia_targets_ch
        def dia_prior_ch
        def dia_rt_window
        if (has_reference_list) {
            STRIP_REFERENCE_DECOYS(reference_list)
            dia_targets_ch = STRIP_REFERENCE_DECOYS.out.reference_list
            dia_prior_ch   = dia_prior
            dia_rt_window  = 0
        } else {
            SELECT_DIATHEM_TARGETS(enc_targets_ch, 'encyclopedia', library, dia_prior,
                                   file("${projectDir}/assets/unimod_delta_masses.tsv"))
            dia_targets_ch = SELECT_DIATHEM_TARGETS.out.targets
            dia_prior_ch   = SELECT_DIATHEM_TARGETS.out.prior.ifEmpty(no_file)
            dia_rt_window  = SELECT_DIATHEM_TARGETS.out.rt_window.map { f -> f.text.trim() }
        }
        DIATHEM_QUANT('DIA', mzmls_ch, dia_targets_ch, dia_prior_ch, sample_map, dia_rt_window)
        diathem_tsv_ch = DIATHEM_QUANT.out.quant.map { _mode, tsv -> tsv }
    } else {
        diathem_tsv_ch = Channel.value(no_file)
    }

    def standard_peak_ch = Channel.value(no_file)
    if (sil && has_reference_list) {
        QUANTIFY_STANDARD_PEAK(samples.map { meta, f, _ml -> tuple(meta, f) }.join(SEARCH_DIA.out.targets)
                                   .combine(tol_ch).map { meta, f, ids, t -> tuple(meta, f, ids, args_of(t, meta.id)) },
                               library, reference_list, heavy_label,
                               'enc_dia')
        standard_peak_ch = QUANTIFY_STANDARD_PEAK.out.rows.collect()
    }
    def coelution_ch = Channel.value(no_file)
    if (sil && has_reference_list) {
        DERIVE_PAIR_MIN_COELUTION(standard_peak_ch, 'enc_dia')
        coelution_ch = DERIVE_PAIR_MIN_COELUTION.out.min_coelution
        if (params.primary_abundance == 'standard_peak') {
            coelution_ch.subscribe { log.info Warnings.minCoelution(it, 'enc_dia') }
        }
    }

    MERGE_QUANT_ENCYCLOPEDIA_DIA(peptides_ch, elib_ch,
                        enc_targets_ch, diathem_tsv_ch, standard_peak_ch, sample_map,
                        heavy_label, params.primary_abundance)


    def tolerance_ch = Channel.value(no_file)
    if (sil) {
        DERIVE_PAIR_RT_TOLERANCE(MERGE_QUANT_ENCYCLOPEDIA_DIA.out.long, 'enc_dia')
        if (params.primary_abundance != 'standard_peak') {
            DERIVE_PAIR_RT_TOLERANCE.out.tolerance.subscribe { log.info Warnings.tolerance(it, 'enc_dia') }
        }
        tolerance_ch = DERIVE_PAIR_RT_TOLERANCE.out.tolerance
    }
    FINALIZE_QUANT(MERGE_QUANT_ENCYCLOPEDIA_DIA.out.long, reference_list, fasta, heavy_label, tolerance_ch, coelution_ch, 'enc_dia')
    warn_ch = warn_ch.mix(FINALIZE_QUANT.out.warnings.flatMap { Warnings.inFile(it, 'enc_dia') })

    def do_abs = sil && standard_amounts.name != 'NO_FILE'
    if (do_abs) {
        QUANTIFY_HEAVY_LIGHT(MERGE_QUANT_ENCYCLOPEDIA_DIA.out.long, standard_amounts, fasta,
                          heavy_label, tolerance_ch, coelution_ch, params.primary_abundance,
                          params.min_consistency, 'enc_dia')
        warn_ch = warn_ch.mix(QUANTIFY_HEAVY_LIGHT.out.warnings
                                    .flatMap { Warnings.inFile(it, 'enc_dia') })
    }

    emit:
    features        = SEARCH_DIA.out.features
    targets         = SEARCH_DIA.out.targets
    decoys          = SEARCH_DIA.out.decoys
    elib            = SEARCH_DIA.out.elib
    base_long       = MERGE_QUANT_ENCYCLOPEDIA_DIA.out.long
    base_per_sample = MERGE_QUANT_ENCYCLOPEDIA_DIA.out.per_sample
    warnings        = warn_ch
}
