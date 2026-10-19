include { OPENSWATH_WORKFLOW          } from '../../modules/local/openms/openswath_workflow'
include { SPLIT_OPENSWATH_FEATURES    } from '../../modules/local/context_ms/split_openswath_features'
include { AUGMENT_CONTEXT_FEATURES    } from '../../modules/local/diathem/augment_context_features'
include { CONTEXT_MS_RUN              } from '../../modules/local/context_ms/context_ms_run'
include { PYPROPHET_SCORE_SINGLE      } from '../../modules/local/pyprophet/pyprophet_score_single'
include { PYPROPHET_TO_PSM            } from '../../modules/local/pyprophet/pyprophet_to_psm'
include { ASSERT_SEARCH_PRODUCTIVE    } from '../../modules/local/context_ms/assert_search_productive'
include { EXTRACT_OPENSWATH_INTENSITY } from '../../modules/local/openms/extract_openswath_intensity'
include { DIATHEM_QUANT               } from '../../modules/local/diathem/diathem_quant'
include { DLIB_TO_REFERENCE_LIST      } from '../../modules/local/encyclopedia/dlib_to_reference_list'
include { STRIP_REFERENCE_DECOYS      } from '../../modules/local/library/strip_reference_decoys'
include { MERGE_QUANT_OPENSWATH_PRM   } from '../../modules/local/quant/merge_quant_openswath_prm'
include { DERIVE_PAIR_RT_TOLERANCE    } from '../../modules/local/quant/derive_pair_rt_tolerance'
include { FINALIZE_QUANT              } from '../../modules/local/quant/finalize_quant'
include { QUANTIFY_HEAVY_LIGHT           } from '../../modules/local/quant/quantify_heavy_light'

workflow OPENSWATH_PRM {
    take:
    samples
    library
    irt_library
    full_window_flag
    spectral_library
    diathem_library
    standard_amounts
    heavy_label
    targets_prm
    no_file
    diathem_ok
    fasta

    main:
    sample_map = samples
        .map { meta, f, _ml -> "${f.baseName}\t${meta.id}" }
        .collectFile(name: 'sample_map.osw_prm.tsv', newLine: true, sort: true)
        .first()

    OPENSWATH_WORKFLOW(samples, library, irt_library, full_window_flag, true)

    def diathem_tsv_ch
    if (diathem_ok) {
        mzmls_ch = samples
            .filter { _m, f, _ml -> f.getName().toLowerCase().endsWith('.mzml') }
            .map    { _m, f, _ml -> f }
            .collect(sort: true)
        def dia_prior = (diathem_library.name != 'NO_FILE') ? diathem_library : spectral_library
        def dia_targets
        if (params.diathem_targets_source == 'library') {
            DLIB_TO_REFERENCE_LIST(spectral_library, 'dlib')
            dia_targets = DLIB_TO_REFERENCE_LIST.out.reference_list
        } else {
            STRIP_REFERENCE_DECOYS(targets_prm)
            dia_targets = STRIP_REFERENCE_DECOYS.out.reference_list
        }
        DIATHEM_QUANT('PRM', mzmls_ch, dia_targets, dia_prior, sample_map, 0)
        diathem_tsv_ch = DIATHEM_QUANT.out.quant.map { _mode, tsv -> tsv }
    } else {
        diathem_tsv_ch = Channel.value(no_file)
    }

    def peptide_ch
    if ((params.openswath_prm_fdr ?: 'context') == 'pyprophet') {
        PYPROPHET_SCORE_SINGLE(OPENSWATH_WORKFLOW.out.osw.map { m, _ml, o -> tuple(m, o) })
        pp_in = samples.map { m, _f, ml -> tuple(m, ml) }
            .join(PYPROPHET_SCORE_SINGLE.out.tsv)
        PYPROPHET_TO_PSM(pp_in)
        peptide_ch = PYPROPHET_TO_PSM.out.peptide
    } else {
        SPLIT_OPENSWATH_FEATURES(OPENSWATH_WORKFLOW.out.osw)

        def split_ch
        if (params.augment_context_with_diathem && diathem_ok) {
            AUGMENT_CONTEXT_FEATURES(SPLIT_OPENSWATH_FEATURES.out.split, diathem_tsv_ch)
            split_ch = AUGMENT_CONTEXT_FEATURES.out.split
        } else {
            split_ch = SPLIT_OPENSWATH_FEATURES.out.split
        }

        CONTEXT_MS_RUN(split_ch)
        peptide_ch = CONTEXT_MS_RUN.out.peptide
    }

    ASSERT_SEARCH_PRODUCTIVE(peptide_ch, 'openswath')
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

    EXTRACT_OPENSWATH_INTENSITY(OPENSWATH_WORKFLOW.out.osw.join(peptide_ch))

    ids_ch   = peptide_ch.map                       { _m, pep -> pep }.collect(sort: true)
    osw_ints = EXTRACT_OPENSWATH_INTENSITY.out.intensity.map { _m, t -> t }.collect(sort: true)
    MERGE_QUANT_OPENSWATH_PRM(ids_ch, osw_ints, diathem_tsv_ch, sample_map,
                        heavy_label, params.primary_abundance)


    def tolerance_ch = Channel.value(no_file)
    if ((heavy_label as List)?.size() > 0) {
        DERIVE_PAIR_RT_TOLERANCE(MERGE_QUANT_OPENSWATH_PRM.out.long, 'osw_prm')
        DERIVE_PAIR_RT_TOLERANCE.out.tolerance.subscribe { log.info Warnings.tolerance(it, 'osw_prm', params.primary_abundance, params.pair_min_coelution) }
        tolerance_ch = DERIVE_PAIR_RT_TOLERANCE.out.tolerance
    }
    FINALIZE_QUANT(MERGE_QUANT_OPENSWATH_PRM.out.long, targets_prm, fasta, heavy_label, tolerance_ch, 'osw_prm')
    def warn_ch = FINALIZE_QUANT.out.warnings
        .flatMap { Warnings.inFile(it, 'osw_prm') }

    def do_abs = (heavy_label as List)?.size() > 0 &&
                 standard_amounts.name != 'NO_FILE'
    if (do_abs) {
        QUANTIFY_HEAVY_LIGHT(MERGE_QUANT_OPENSWATH_PRM.out.long, standard_amounts, fasta,
                          heavy_label, tolerance_ch, params.primary_abundance,
                          params.min_consistency, 'osw_prm')
        warn_ch = warn_ch.mix(QUANTIFY_HEAVY_LIGHT.out.warnings
                                    .flatMap { Warnings.inFile(it, 'osw_prm') })
    }

    emit:
    base_long       = MERGE_QUANT_OPENSWATH_PRM.out.long
    base_per_sample = MERGE_QUANT_OPENSWATH_PRM.out.per_sample
    warnings        = warn_ch
}
