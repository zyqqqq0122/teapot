include { ASSERT_OPENSWATH_FEASIBLE   } from '../../modules/local/openms/assert_openswath_feasible'
include { OPENSWATH_WORKFLOW          } from '../../modules/local/openms/openswath_workflow'
include { SPLIT_OPENSWATH_FEATURES    } from '../../modules/local/context_ms/split_openswath_features'
include { AUGMENT_CONTEXT_FEATURES    } from '../../modules/local/diathem/augment_context_features'
include { CONTEXT_MS_RUN              } from '../../modules/local/context_ms/context_ms_run'
include { PYPROPHET_SCORE_SINGLE      } from '../../modules/local/pyprophet/pyprophet_score_single'
include { PYPROPHET_TO_PSM            } from '../../modules/local/pyprophet/pyprophet_to_psm'
include { ASSERT_SEARCH_PRODUCTIVE    } from '../../modules/local/context_ms/assert_search_productive'
include { EXTRACT_OPENSWATH_INTENSITY } from '../../modules/local/openms/extract_openswath_intensity'
include { DIATHEM_QUANT               } from '../../modules/local/diathem/diathem_quant'
include { STRIP_REFERENCE_DECOYS      } from '../../modules/local/library/strip_reference_decoys'
include { MERGE_QUANT_OPENSWATH_PRM   } from '../../modules/local/quant/merge_quant_openswath_prm'
include { FINALIZE_QUANT              } from '../../modules/local/quant/finalize_quant'
include { QUANTIFY_HEAVY_LIGHT           } from '../../modules/local/quant/quantify_heavy_light'
include { CONTEXT_SEARCH        as OSWPRM_CONTEXT_SEARCH    } from '../../modules/local/encyclopedia/context_search'
include { ENCYCLOPEDIA_LIBEXPORT as OPENSWATH_PRM_LIBEXPORT } from '../../modules/local/encyclopedia/encyclopedia_libexport'

workflow OPENSWATH_PRM {
    take:
    samples
    library
    irt_library
    full_window_flag
    encyclopedia_library
    fasta
    diathem_library
    standard_amounts
    heavy_label
    sample_map
    targets_prm
    no_file
    diathem_ok

    main:
    def nested_ch
    if (!params.skip_openswath_feasibility) {
        ASSERT_OPENSWATH_FEASIBLE(samples.map { meta, mzml, _rl -> tuple(meta, mzml) })
        nested_ch = ASSERT_OPENSWATH_FEASIBLE.out.nested
            .map { _m, f -> f }
            .first()
            .ifEmpty(no_file)
    } else {
        nested_ch = Channel.value(no_file)
    }

    OPENSWATH_WORKFLOW(samples, library, irt_library, full_window_flag, nested_ch)
    EXTRACT_OPENSWATH_INTENSITY(OPENSWATH_WORKFLOW.out.osw)

    def enc_elib_ch
    if (params.openswath_prm_encyclopedia_quant) {
        OSWPRM_CONTEXT_SEARCH(samples, encyclopedia_library, fasta)
        pooled_enc = OSWPRM_CONTEXT_SEARCH.out.search_products
            .map { _m, files -> files }
            .flatten()
            .collect(sort: true)
        OPENSWATH_PRM_LIBEXPORT(pooled_enc, 'openswath_prm_encyclopedia', encyclopedia_library, fasta, params.context_search_args)
        enc_elib_ch = OPENSWATH_PRM_LIBEXPORT.out.elib
    } else {
        enc_elib_ch = Channel.value(no_file)
    }

    def diathem_tsv_ch
    if (diathem_ok) {
        mzmls_ch = samples
            .filter { _m, f, _ml -> f.getName().toLowerCase() ==~ /.*\.(mzml|mzxml)$/ }
            .map    { _m, f, _ml -> f }
            .collect(sort: true)
        def dia_prior = (diathem_library.name != 'NO_FILE') ? diathem_library : library
        def dia_targets
        if (params.diathem_targets_source == 'library') {
            dia_targets = library
        } else {
            STRIP_REFERENCE_DECOYS(targets_prm)
            dia_targets = STRIP_REFERENCE_DECOYS.out.reference_list
        }
        DIATHEM_QUANT('PRM', mzmls_ch, dia_targets, dia_prior, sample_map)
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

    ids_ch   = peptide_ch.map                       { _m, pep -> pep }.collect(sort: true)
    osw_ints = EXTRACT_OPENSWATH_INTENSITY.out.intensity.map { _m, t -> t }.collect(sort: true)
    MERGE_QUANT_OPENSWATH_PRM(ids_ch, osw_ints, diathem_tsv_ch, sample_map,
                        heavy_label, params.primary_abundance, enc_elib_ch)


    FINALIZE_QUANT(MERGE_QUANT_OPENSWATH_PRM.out.long, targets_prm, heavy_label, 'osw_prm')

    def do_abs = (heavy_label as List)?.size() > 0 &&
                 standard_amounts.name != 'NO_FILE'
    if (do_abs) {
        QUANTIFY_HEAVY_LIGHT(MERGE_QUANT_OPENSWATH_PRM.out.long, standard_amounts,
                          heavy_label, params.primary_abundance,
                          params.min_consistency)
    }

    emit:
    base_long       = MERGE_QUANT_OPENSWATH_PRM.out.long
    base_per_sample = MERGE_QUANT_OPENSWATH_PRM.out.per_sample
}
