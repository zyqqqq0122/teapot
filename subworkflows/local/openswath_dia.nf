include { OPENSWATH_WORKFLOW        } from '../../modules/local/openms/openswath_workflow'
include { PYPROPHET_SUBSAMPLE       } from '../../modules/local/pyprophet/pyprophet_subsample'
include { PYPROPHET_MERGE           } from '../../modules/local/pyprophet/pyprophet_merge'
include { PYPROPHET_LEARN           } from '../../modules/local/pyprophet/pyprophet_learn'
include { PYPROPHET_APPLY           } from '../../modules/local/pyprophet/pyprophet_apply'
include { PYPROPHET_EXPORT          } from '../../modules/local/pyprophet/pyprophet_export'
include { DIATHEM_QUANT             } from '../../modules/local/diathem/diathem_quant'
include { MERGE_QUANT_OPENSWATH_DIA } from '../../modules/local/quant/merge_quant_openswath_dia'
include { FINALIZE_QUANT            } from '../../modules/local/quant/finalize_quant'
include { QUANTIFY_HEAVY_LIGHT         } from '../../modules/local/quant/quantify_heavy_light'
include { TRIC_FEATURE_ALIGNMENT    } from '../../modules/local/tric/tric_feature_alignment'

workflow OPENSWATH_DIA {
    take:
    samples
    library
    irt_library
    full_window_flag
    fasta
    diathem_library
    standard_amounts
    heavy_label
    sample_map
    no_file
    diathem_ok

    main:
    def nested_ch = Channel.value(no_file)

    OPENSWATH_WORKFLOW(samples, library, irt_library, full_window_flag, nested_ch)

    n_runs_ch = OPENSWATH_WORKFLOW.out.osw.count()
    PYPROPHET_SUBSAMPLE(OPENSWATH_WORKFLOW.out.osw, n_runs_ch)

    subsampled = PYPROPHET_SUBSAMPLE.out.osws.map { _m, o -> o }.collect(sort: true)
    PYPROPHET_MERGE(subsampled, library)
    PYPROPHET_LEARN(PYPROPHET_MERGE.out.merged)

    PYPROPHET_APPLY(OPENSWATH_WORKFLOW.out.osw, PYPROPHET_LEARN.out.model)
    PYPROPHET_EXPORT(PYPROPHET_APPLY.out.scored)

    def diathem_tsv_ch
    if (diathem_ok) {
        mzmls_ch = samples
            .filter { _m, f, _ml -> f.getName().toLowerCase() ==~ /.*\.(mzml|mzxml)$/ }
            .map    { _m, f, _ml -> f }
            .collect(sort: true)
        def dia_prior = (diathem_library.name != 'NO_FILE') ? diathem_library : library
        DIATHEM_QUANT('DIA', mzmls_ch, library, dia_prior, sample_map)
        diathem_tsv_ch = DIATHEM_QUANT.out.quant.map { _mode, tsv -> tsv }
    } else {
        diathem_tsv_ch = Channel.value(no_file)
    }

    pp_tsvs = PYPROPHET_EXPORT.out.tsv.map { _m, _ml, tsv -> tsv }.collect(sort: true)

    def tric_tsv_ch
    if (params.run_tric) {
        TRIC_FEATURE_ALIGNMENT(pp_tsvs)
        tric_tsv_ch = TRIC_FEATURE_ALIGNMENT.out.long
    } else {
        tric_tsv_ch = Channel.value(no_file)
    }

    MERGE_QUANT_OPENSWATH_DIA(pp_tsvs, diathem_tsv_ch, tric_tsv_ch, sample_map,
                        heavy_label, params.primary_abundance)


    FINALIZE_QUANT(MERGE_QUANT_OPENSWATH_DIA.out.long, no_file, heavy_label, 'osw_dia')

    def do_abs = (heavy_label as List)?.size() > 0 &&
                 standard_amounts.name != 'NO_FILE'
    if (do_abs) {
        QUANTIFY_HEAVY_LIGHT(MERGE_QUANT_OPENSWATH_DIA.out.long, standard_amounts,
                          heavy_label, params.primary_abundance,
                          params.min_consistency, 'osw_dia')
    }

    emit:
    base_long       = MERGE_QUANT_OPENSWATH_DIA.out.long
    base_per_sample = MERGE_QUANT_OPENSWATH_DIA.out.per_sample
}
