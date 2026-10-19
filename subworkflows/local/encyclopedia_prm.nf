include { CONTEXT_SEARCH               } from '../../modules/local/encyclopedia/context_search'
include { ASSERT_SEARCH_PRODUCTIVE     } from '../../modules/local/context_ms/assert_search_productive'
include { ENCYCLOPEDIA_LIBEXPORT       } from '../../modules/local/encyclopedia/encyclopedia_libexport'
include { ENCYCLOPEDIA_LIBEXPORT_PER_SAMPLE         } from '../../modules/local/encyclopedia/encyclopedia_libexport_per_sample'
include { ASSERT_QUANT_PRODUCTIVE      } from '../../modules/local/quant/assert_quant_productive'
include { DIATHEM_QUANT                } from '../../modules/local/diathem/diathem_quant'
include { STRIP_REFERENCE_DECOYS      } from '../../modules/local/library/strip_reference_decoys'
include { MERGE_QUANT_ENCYCLOPEDIA_PRM } from '../../modules/local/quant/merge_quant_encyclopedia_prm'
include { DERIVE_PAIR_RT_TOLERANCE     } from '../../modules/local/quant/derive_pair_rt_tolerance'
include { FINALIZE_QUANT               } from '../../modules/local/quant/finalize_quant'
include { QUANTIFY_HEAVY_LIGHT            } from '../../modules/local/quant/quantify_heavy_light'

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

    main:
    sample_map = samples
        .map { meta, f, _ml -> "${f.baseName}\t${meta.id}" }
        .collectFile(name: 'sample_map.enc_prm.tsv', newLine: true, sort: true)
        .first()

    CONTEXT_SEARCH(samples, elib, fasta)

    ASSERT_SEARCH_PRODUCTIVE(CONTEXT_SEARCH.out.ref_targets, 'encyclopedia')
    ASSERT_SEARCH_PRODUCTIVE.out.qc
        .map { _m, qc -> def l = qc.readLines(); l.size() > 1 ? (l[1].tokenize('\t')[4] as int) : 1 }
        .collect()
        .subscribe { n ->
            if (n.sum() == 0) {
                error "ASSERT_SEARCH_PRODUCTIVE: no reference peptide at q <= ${params.fdr} in any " +
                      "of the ${n.size()} runs. Check the search tolerances for your instrument " +
                      "first; qc/<sample>/ holds each run's report."
            }
        }

    ENCYCLOPEDIA_LIBEXPORT_PER_SAMPLE(CONTEXT_SEARCH.out.search_products, elib, fasta,
                         params.context_search_args)

    pooled_prm = CONTEXT_SEARCH.out.search_products
        .map { _m, files -> files }
        .flatten()
        .mix(ENCYCLOPEDIA_LIBEXPORT_PER_SAMPLE.out.elib.map { _m, e -> e }.flatten())
        .collect(sort: true)
    ENCYCLOPEDIA_LIBEXPORT(pooled_prm, 'analysis', elib, fasta, params.context_search_args, true)

    ASSERT_QUANT_PRODUCTIVE(ENCYCLOPEDIA_LIBEXPORT.out.elib, 'encyclopedia_prm')

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

    MERGE_QUANT_ENCYCLOPEDIA_PRM(ref_tgts, ref_feats, diathem_tsv_ch, sample_map,
                        heavy_label, params.primary_abundance,
                        ENCYCLOPEDIA_LIBEXPORT.out.elib)


    def tolerance_ch = Channel.value(no_file)
    if ((heavy_label as List)?.size() > 0) {
        DERIVE_PAIR_RT_TOLERANCE(MERGE_QUANT_ENCYCLOPEDIA_PRM.out.long, 'enc_prm')
        DERIVE_PAIR_RT_TOLERANCE.out.tolerance.subscribe { log.info Warnings.tolerance(it, 'enc_prm') }
        tolerance_ch = DERIVE_PAIR_RT_TOLERANCE.out.tolerance
    }
    FINALIZE_QUANT(MERGE_QUANT_ENCYCLOPEDIA_PRM.out.long, targets_prm, fasta, heavy_label, tolerance_ch, 'enc_prm')
    def warn_ch = ASSERT_QUANT_PRODUCTIVE.out.warnings
        .mix(FINALIZE_QUANT.out.warnings)
        .flatMap { Warnings.inFile(it, 'enc_prm') }

    def do_abs = (heavy_label as List)?.size() > 0 &&
                 standard_amounts.name != 'NO_FILE'
    if (do_abs) {
        QUANTIFY_HEAVY_LIGHT(MERGE_QUANT_ENCYCLOPEDIA_PRM.out.long, standard_amounts, fasta,
                          heavy_label, tolerance_ch, params.primary_abundance,
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
