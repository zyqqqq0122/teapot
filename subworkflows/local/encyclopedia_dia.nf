include { SEARCH_DIA                   } from '../../modules/local/encyclopedia/search_dia'
include { ENCYCLOPEDIA_LIBEXPORT       } from '../../modules/local/encyclopedia/encyclopedia_libexport'
include { ASSERT_QUANT_PRODUCTIVE      } from '../../modules/local/quant/assert_quant_productive'
include { DIATHEM_QUANT                } from '../../modules/local/diathem/diathem_quant'
include { SELECT_DIATHEM_TARGETS       } from '../../modules/local/diathem/select_diathem_targets'
include { STRIP_REFERENCE_DECOYS      } from '../../modules/local/library/strip_reference_decoys'
include { MERGE_QUANT_ENCYCLOPEDIA_DIA } from '../../modules/local/quant/merge_quant_encyclopedia_dia'
include { FINALIZE_QUANT               } from '../../modules/local/quant/finalize_quant'
include { QUANTIFY_HEAVY_LIGHT            } from '../../modules/local/quant/quantify_heavy_light'

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

    main:
    sample_map = samples
        .map { meta, f, _ml -> "${f.baseName}\t${meta.id}" }
        .collectFile(name: 'sample_map.enc_dia.tsv', newLine: true, sort: true)
        .first()

    dia_in = samples.map { meta, dia, _ml -> tuple(meta, dia) }
    SEARCH_DIA(dia_in, library, fasta)

    pooled = SEARCH_DIA.out.features.map { _m, f -> f }
        .mix(SEARCH_DIA.out.targets.map { _m, t -> t })
        .mix(SEARCH_DIA.out.decoys .map { _m, d -> d })
        .mix(SEARCH_DIA.out.elib   .map { _m, e -> e })
        .mix(samples.map { _m, f, _ml -> f })
        .collect(sort: true)
    ENCYCLOPEDIA_LIBEXPORT(pooled, 'analysis', library, fasta, params.dia_search_args,
                           has_reference_list)

    ASSERT_QUANT_PRODUCTIVE(ENCYCLOPEDIA_LIBEXPORT.out.elib, 'encyclopedia_dia')

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

    peptides_ch = ENCYCLOPEDIA_LIBEXPORT.out.peptides
        .ifEmpty {
            log.warn "ENCYCLOPEDIA_DIA: libexport wrote no analysis.elib.peptides.txt " +
                     "(single search, no cross-run alignment). Abundance comes " +
                     "from the .elib peptidequants table, which is the primary " +
                     "source in either case."
            return file("${projectDir}/assets/NO_FILE")
        }
    MERGE_QUANT_ENCYCLOPEDIA_DIA(peptides_ch, ENCYCLOPEDIA_LIBEXPORT.out.elib,
                        enc_targets_ch, diathem_tsv_ch, sample_map,
                        heavy_label, params.primary_abundance)


    FINALIZE_QUANT(MERGE_QUANT_ENCYCLOPEDIA_DIA.out.long, reference_list, fasta, heavy_label, 'enc_dia')

    def do_abs = (heavy_label as List)?.size() > 0 &&
                 standard_amounts.name != 'NO_FILE'
    if (do_abs) {
        QUANTIFY_HEAVY_LIGHT(MERGE_QUANT_ENCYCLOPEDIA_DIA.out.long, standard_amounts, fasta,
                          heavy_label, params.primary_abundance,
                          params.min_consistency, 'enc_dia')
    }

    emit:
    features        = SEARCH_DIA.out.features
    targets         = SEARCH_DIA.out.targets
    decoys          = SEARCH_DIA.out.decoys
    elib            = SEARCH_DIA.out.elib
    base_long       = MERGE_QUANT_ENCYCLOPEDIA_DIA.out.long
    base_per_sample = MERGE_QUANT_ENCYCLOPEDIA_DIA.out.per_sample
}
