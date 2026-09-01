include { SEARCH_DIA                   } from '../../modules/local/encyclopedia/search_dia'
include { ENCYCLOPEDIA_LIBEXPORT       } from '../../modules/local/encyclopedia/encyclopedia_libexport'
include { ASSERT_QUANT_PRODUCTIVE      } from '../../modules/local/quant/assert_quant_productive'
include { DIATHEM_QUANT                } from '../../modules/local/diathem/diathem_quant'
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
    sample_map
    diathem_targets
    no_file
    diathem_ok

    main:
    dia_in = samples.map { meta, dia, _ml -> tuple(meta, dia) }
    SEARCH_DIA(dia_in, library, fasta)

    pooled = SEARCH_DIA.out.features.map { _m, f -> f }
        .mix(SEARCH_DIA.out.targets.map { _m, t -> t })
        .mix(SEARCH_DIA.out.decoys .map { _m, d -> d })
        .mix(SEARCH_DIA.out.elib   .map { _m, e -> e })
        .mix(samples.map { _m, f, _ml -> f })
        .collect(sort: true)
    ENCYCLOPEDIA_LIBEXPORT(pooled, 'analysis', library, fasta, params.dia_search_args)

    ASSERT_QUANT_PRODUCTIVE(ENCYCLOPEDIA_LIBEXPORT.out.elib, 'encyclopedia_dia')

    def diathem_tsv_ch
    def have_targets = diathem_targets.name != 'NO_FILE'
    if (diathem_ok && have_targets) {
        mzmls_ch = samples
            .filter { _m, f, _ml -> f.getName().toLowerCase() ==~ /.*\.(mzml|mzxml)$/ }
            .map    { _m, f, _ml -> f }
            .collect(sort: true)
        def dia_prior = (diathem_library.name != 'NO_FILE') ? diathem_library : library
        DIATHEM_QUANT('DIA', mzmls_ch, diathem_targets, dia_prior, sample_map)
        diathem_tsv_ch = DIATHEM_QUANT.out.quant.map { _mode, tsv -> tsv }
    } else {
        if (diathem_ok && !have_targets) {
            log.warn "ENCYCLOPEDIA_DIA: skipping DIATHEM_QUANT, no target list. " +
                     "diathem needs a Skyline/EncyclopeDIA assay (.csv/.tsv/.txt) " +
                     "or an OpenSWATH .pqp; a .dlib/.elib is not readable as one. " +
                     "Supply the samplesheet's reference_list column, or --blib so " +
                     "BLIB_TO_REFERENCE_LIST can derive it."
        }
        diathem_tsv_ch = Channel.value(no_file)
    }

    enc_targets_ch = SEARCH_DIA.out.targets.map { _m, t -> t }.collect(sort: true)

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


    FINALIZE_QUANT(MERGE_QUANT_ENCYCLOPEDIA_DIA.out.long, diathem_targets, heavy_label, 'enc_dia')

    def do_abs = (heavy_label as List)?.size() > 0 &&
                 standard_amounts.name != 'NO_FILE'
    if (do_abs) {
        QUANTIFY_HEAVY_LIGHT(MERGE_QUANT_ENCYCLOPEDIA_DIA.out.long, standard_amounts,
                          heavy_label, params.primary_abundance,
                          params.min_consistency)
    }

    emit:
    features        = SEARCH_DIA.out.features
    targets         = SEARCH_DIA.out.targets
    decoys          = SEARCH_DIA.out.decoys
    elib            = SEARCH_DIA.out.elib
    base_long       = MERGE_QUANT_ENCYCLOPEDIA_DIA.out.long
    base_per_sample = MERGE_QUANT_ENCYCLOPEDIA_DIA.out.per_sample
}
