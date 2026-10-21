include { CONTEXT_SEARCH               } from '../../modules/local/encyclopedia/context_search'
include { ALIGN_RETENTION_TIMES        } from '../../modules/local/runs/align_retention_times'
include { ASSERT_SEARCH_PRODUCTIVE     } from '../../modules/local/context_ms/assert_search_productive'
include { ENCYCLOPEDIA_LIBEXPORT       } from '../../modules/local/encyclopedia/encyclopedia_libexport'
include { ENCYCLOPEDIA_LIBEXPORT_PER_SAMPLE         } from '../../modules/local/encyclopedia/encyclopedia_libexport_per_sample'
include { ASSERT_QUANT_PRODUCTIVE      } from '../../modules/local/quant/assert_quant_productive'
include { DIATHEM_QUANT                } from '../../modules/local/diathem/diathem_quant'
include { STRIP_REFERENCE_DECOYS      } from '../../modules/local/library/strip_reference_decoys'
include { MERGE_QUANT_ENCYCLOPEDIA_PRM } from '../../modules/local/quant/merge_quant_encyclopedia_prm'
include { DERIVE_PAIR_RT_TOLERANCE     } from '../../modules/local/quant/derive_pair_rt_tolerance'
include { DERIVE_PAIR_MIN_COELUTION    } from '../../modules/local/quant/derive_pair_min_coelution'
include { DERIVE_SEARCH_TOLERANCE      } from '../../modules/local/runs/derive_search_tolerance'
include { FINALIZE_QUANT               } from '../../modules/local/quant/finalize_quant'
include { QUANTIFY_HEAVY_LIGHT            } from '../../modules/local/quant/quantify_heavy_light'
include { QUANTIFY_STANDARD_PEAK       } from '../../modules/local/quant/quantify_standard_peak'

workflow ENCYCLOPEDIA_PRM {
    take:
    samples
    elib
    fasta
    diathem_library
    standard_amounts
    heavy_label
    targets_prm
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
        .collectFile(name: 'sample_map.enc_prm.tsv', newLine: true, sort: true)
        .first()

    def sil = (heavy_label as List)?.size() > 0
    def warn_ch = Channel.empty()

    def runs_in
    if (rt_alignment && sil) {
        def ids = samples.map { meta, _f, _ml -> meta.id }
            .collectFile(name: 'samples.enc_prm.txt', newLine: true, sort: true)
        def mine = samples.map { meta, _f, _ml -> [meta.id] }
            .join(located.map { meta, t -> [meta.id, t] })
            .map { _id, t -> t }
        ALIGN_RETENTION_TIMES(ids, mine.collect().ifEmpty([]), acquired.collect().ifEmpty([]),
                              elib, targets_prm, heavy_label, 'enc_prm')
        warn_ch = warn_ch.mix(ALIGN_RETENTION_TIMES.out.warnings.flatMap { Warnings.inFile(it, 'enc_prm') })
        def windows = ALIGN_RETENTION_TIMES.out.windows.flatten()
            .map { w -> [w.name - '.apex_windows.tsv', w] }
        runs_in = samples.map { meta, f, rl -> [meta.id, meta, f, rl] }
            .join(windows, remainder: true)
            .map { _id, meta, f, rl, w -> tuple(meta, f, rl, w ?: no_file) }
    } else {
        runs_in = samples.map { meta, f, rl -> tuple(meta, f, rl, no_file) }
    }
    DERIVE_SEARCH_TOLERANCE(fragment_errors.map { meta, f -> tuple(meta.id, f) }
                                .join(samples.map { meta, _f, _ml -> tuple(meta.id, true) })
                                .map { _id, f, _s -> f }.collect().ifEmpty([no_file]),
                            Args.fragmentTolerance(params.context_search_args), 'enc_prm')
    def tol_ch = DERIVE_SEARCH_TOLERANCE.out.tolerance.combine(DERIVE_SEARCH_TOLERANCE.out.runs)
        .map { t, r -> Args.readTolerance(t, r) }
        .first()
    def args_of = { Map t, String sid -> t.runs.contains(sid) ? Args.withBatchTolerance(params.context_search_args, t) : Args.flat(params.context_search_args) }
    def batch_args = tol_ch.map { t -> Args.withBatchTolerance(params.context_search_args, t) }
    tol_ch.subscribe { t -> log.info Warnings.fragmentTolerance(t, 'enc_prm', params.context_search_args) }
    CONTEXT_SEARCH(runs_in.combine(tol_ch).map { meta, f, rl, w, t -> tuple(meta, f, rl, w, args_of(t, meta.id)) }, elib, fasta)

    ASSERT_SEARCH_PRODUCTIVE(CONTEXT_SEARCH.out.ref_targets, 'encyclopedia')
    ASSERT_SEARCH_PRODUCTIVE.out.qc
        .map { _m, qc -> def l = qc.readLines(); l.size() > 1 ? (l[1].tokenize('\t')[4] as int) : 1 }
        .collect()
        .ifEmpty([])
        .subscribe { n ->
            if (!n) {
                error "ENCYCLOPEDIA_PRM: the search failed for every run; the errors above say why."
            }
            if (n.sum() == 0) {
                error "ASSERT_SEARCH_PRODUCTIVE: no reference peptide at q <= ${params.fdr} in any " +
                      "of the ${n.size()} runs. Check the search tolerances for your instrument " +
                      "first; qc/<sample>/ holds each run's report."
            }
        }

    def quant_elib_ch = Channel.value(no_file)
    if (enc_quant) {
        ENCYCLOPEDIA_LIBEXPORT_PER_SAMPLE(CONTEXT_SEARCH.out.search_products, elib, fasta,
                             batch_args)

        pooled_prm = CONTEXT_SEARCH.out.search_products
            .map { _m, files -> files }
            .flatten()
            .mix(ENCYCLOPEDIA_LIBEXPORT_PER_SAMPLE.out.elib.map { _m, e -> e }.flatten())
            .collect(sort: true)
        ENCYCLOPEDIA_LIBEXPORT(pooled_prm, 'analysis', elib, fasta, batch_args, true)

        ASSERT_QUANT_PRODUCTIVE(ENCYCLOPEDIA_LIBEXPORT.out.elib, 'encyclopedia_prm')
        warn_ch = warn_ch.mix(ASSERT_QUANT_PRODUCTIVE.out.warnings.flatMap { Warnings.inFile(it, 'enc_prm') })
        quant_elib_ch = ENCYCLOPEDIA_LIBEXPORT.out.elib
    }

    def diathem_tsv_ch
    if (diathem_ok) {
        mzmls_ch = samples
            .filter { _m, f, _ml -> f.getName().toLowerCase().endsWith('.mzml') }
            .map    { _m, f, _ml -> f }
            .collect(sort: true)
        def dia_prior = (diathem_library.name != 'NO_FILE') ? diathem_library : elib
        STRIP_REFERENCE_DECOYS(targets_prm)
        DIATHEM_QUANT('PRM', mzmls_ch, STRIP_REFERENCE_DECOYS.out.reference_list,
                      dia_prior, sample_map, 0)
        diathem_tsv_ch = DIATHEM_QUANT.out.quant.map { _mode, tsv -> tsv }
    } else {
        diathem_tsv_ch = Channel.value(no_file)
    }

    ref_tgts  = CONTEXT_SEARCH.out.ref_targets.map  { _m, t -> t }.collect(sort: true)
    ref_feats = CONTEXT_SEARCH.out.ref_features.map { _m, t -> t }.collect(sort: true)

    def standard_peak_ch = Channel.value(no_file)
    if (sil) {
        QUANTIFY_STANDARD_PEAK(samples.map { meta, f, _rl -> tuple(meta, f) }.join(CONTEXT_SEARCH.out.ref_targets)
                                   .combine(tol_ch).map { meta, f, ids, t -> tuple(meta, f, ids, args_of(t, meta.id)) },
                               elib, targets_prm, heavy_label,
                               'enc_prm')
        standard_peak_ch = QUANTIFY_STANDARD_PEAK.out.rows.collect()
    }
    def coelution_ch = Channel.value(no_file)
    if (sil) {
        DERIVE_PAIR_MIN_COELUTION(standard_peak_ch, 'enc_prm')
        coelution_ch = DERIVE_PAIR_MIN_COELUTION.out.min_coelution
        if (params.primary_abundance == 'standard_peak') {
            coelution_ch.subscribe { log.info Warnings.minCoelution(it, 'enc_prm') }
        }
    }

    MERGE_QUANT_ENCYCLOPEDIA_PRM(ref_tgts, ref_feats, diathem_tsv_ch, standard_peak_ch, sample_map,
                        heavy_label, params.primary_abundance,
                        quant_elib_ch)


    def tolerance_ch = Channel.value(no_file)
    if (sil) {
        DERIVE_PAIR_RT_TOLERANCE(MERGE_QUANT_ENCYCLOPEDIA_PRM.out.long, 'enc_prm')
        if (params.primary_abundance != 'standard_peak') {
            DERIVE_PAIR_RT_TOLERANCE.out.tolerance.subscribe { log.info Warnings.tolerance(it, 'enc_prm') }
        }
        tolerance_ch = DERIVE_PAIR_RT_TOLERANCE.out.tolerance
    }
    FINALIZE_QUANT(MERGE_QUANT_ENCYCLOPEDIA_PRM.out.long, targets_prm, fasta, heavy_label, tolerance_ch, coelution_ch,
                   sample_map, 'enc_prm')
    warn_ch = warn_ch.mix(FINALIZE_QUANT.out.warnings.flatMap { Warnings.inFile(it, 'enc_prm') })

    def do_abs = sil && standard_amounts.name != 'NO_FILE'
    if (do_abs) {
        QUANTIFY_HEAVY_LIGHT(MERGE_QUANT_ENCYCLOPEDIA_PRM.out.long, standard_amounts, fasta,
                          heavy_label, tolerance_ch, coelution_ch, params.primary_abundance,
                          params.min_consistency, 'enc_prm')
        warn_ch = warn_ch.mix(QUANTIFY_HEAVY_LIGHT.out.warnings
                                    .flatMap { Warnings.inFile(it, 'enc_prm') })
    }

    emit:
    reference_targets  = CONTEXT_SEARCH.out.ref_targets
    reference_decoys   = CONTEXT_SEARCH.out.ref_decoys
    background_targets = CONTEXT_SEARCH.out.bg_targets
    background_decoys  = CONTEXT_SEARCH.out.bg_decoys
    plots              = CONTEXT_SEARCH.out.plots
    base_long          = MERGE_QUANT_ENCYCLOPEDIA_PRM.out.long
    base_per_sample    = MERGE_QUANT_ENCYCLOPEDIA_PRM.out.per_sample
    warnings           = warn_ch
}
