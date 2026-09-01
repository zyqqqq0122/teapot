include { KOINA_PREDICT                           } from '../../modules/local/encyclopedia/koina_predict'
include { BLIB_TO_DLIB                            } from '../../modules/local/encyclopedia/blib_to_dlib'
include { BLIB_TO_REFERENCE_LIST                  } from '../../modules/local/encyclopedia/blib_to_reference_list'
include { ADD_DECOYS_TO_REFERENCE_LIST               } from '../../modules/local/encyclopedia/add_decoys_to_reference_list'
include { COMPLETE_CHANNELS                        } from '../../modules/local/library/complete_channels'
include { EMIT_IRT_LIBRARY                        } from '../../modules/local/library/emit_irt_library'
include { HARMONIZE_RT as HARMONIZE_RT_TARGET     } from '../../modules/local/library/harmonize_rt'
include { HARMONIZE_RT as HARMONIZE_RT_BACKGROUND } from '../../modules/local/library/harmonize_rt'
include { MERGE_LIBRARIES                           } from '../../modules/local/library/merge_libraries'
include { ASSERT_RT_COHERENT                      } from '../../modules/local/library/assert_rt_coherent'
include { ASSERT_BACKGROUND_ADEQUATE              } from '../../modules/local/library/assert_background_adequate'
include { MERGE_GPF                               } from '../../modules/local/encyclopedia/merge_gpf'
include { SEARCH_GPF                              } from '../../modules/local/encyclopedia/search_gpf'
include { ENCYCLOPEDIA_LIBEXPORT                  } from '../../modules/local/encyclopedia/encyclopedia_libexport'
include { PREPARE_BACKGROUND                      } from './prepare_background'

workflow PREPARE_LIBRARY_ENCYCLOPEDIA {
    take:
    fasta
    heavy_label
    background_source
    background_library
    background_fasta
    background_min_targets
    no_file

    main:
    def gpf_is_background = (background_source == 'gpf')

    def need_dlib   = (params.library_sheet && !gpf_is_background) || !params.elib
    def run_channel = (heavy_label as List)?.size() > 0
    def has_bg      = background_source && background_source != 'none'
    def no_file_ch  = Channel.value(no_file)

    def dlib_ch = null
    def reference_list_ch = no_file_ch

    if (need_dlib) {
        if (params.dlib) {
            dlib_ch = Channel.value(file(params.dlib))
        } else if (params.blib) {
            def blib_ch = Channel.value(file(params.blib))
            def source_blib_ch
            if (run_channel) {
                COMPLETE_CHANNELS(blib_ch, heavy_label)
                source_blib_ch = COMPLETE_CHANNELS.out.library
            } else {
                source_blib_ch = blib_ch
            }
            BLIB_TO_DLIB(source_blib_ch, fasta)

            BLIB_TO_REFERENCE_LIST(source_blib_ch)

            if (params.add_decoys_to_reference_list) {
                ADD_DECOYS_TO_REFERENCE_LIST(BLIB_TO_REFERENCE_LIST.out.reference_list,
                                     Channel.value(file(params.reference_list_decoy_jar)))
                reference_list_ch = ADD_DECOYS_TO_REFERENCE_LIST.out.reference_list
            } else {
                reference_list_ch = BLIB_TO_REFERENCE_LIST.out.reference_list
            }

            if (has_bg) {
                def anchors_ch
                if (params.irt_library && params.irt_anchors) {
                    anchors_ch = Channel.value(file(params.irt_anchors))
                } else {
                    EMIT_IRT_LIBRARY(blib_ch, heavy_label)
                    anchors_ch = EMIT_IRT_LIBRARY.out.anchors
                }
                HARMONIZE_RT_TARGET(BLIB_TO_DLIB.out.dlib, anchors_ch,
                                    'target', false, false)

                PREPARE_BACKGROUND(fasta, background_source,
                                   background_library, background_fasta)
                HARMONIZE_RT_BACKGROUND(PREPARE_BACKGROUND.out.dlib,
                                        anchors_ch, 'background',
                                        false, false)

                MERGE_LIBRARIES(HARMONIZE_RT_TARGET.out.library,
                              HARMONIZE_RT_BACKGROUND.out.library)
                dlib_ch = MERGE_LIBRARIES.out.dlib

                ASSERT_RT_COHERENT(dlib_ch)
                ASSERT_BACKGROUND_ADEQUATE(dlib_ch, reference_list_ch,
                                           background_min_targets)
            } else {
                if (run_channel) {
                    log.warn "PREPARE_LIBRARY_ENCYCLOPEDIA: --blib + heavy_label without a " +
                             "background library will produce reference_list == " +
                             "library and CONTEXT_SEARCH will find zero background " +
                             "targets. Set background_source to 'koina' | 'dlib' | " +
                             "'blib' | 'elib' | 'gpf'."
                }
                dlib_ch = BLIB_TO_DLIB.out.dlib
            }
        } else if (params.use_koina) {
            KOINA_PREDICT(fasta)
            def source_dlib_ch
            if (run_channel) {
                COMPLETE_CHANNELS(KOINA_PREDICT.out.dlib, heavy_label)
                source_dlib_ch = COMPLETE_CHANNELS.out.library
            } else {
                source_dlib_ch = KOINA_PREDICT.out.dlib
            }
            dlib_ch = source_dlib_ch
        } else {
            error "No spectral library available. Provide --dlib, --blib, --use_koina, " +
                  "or supply --elib to skip .dlib entirely when not building from GPF."
        }
    }

    def library_ch
    if (params.elib && (!params.library_sheet || gpf_is_background)) {
        library_ch = Channel.value(file(params.elib))

    } else if (params.library_sheet && !gpf_is_background) {
        rows = Channel
            .fromPath(params.library_sheet)
            .splitCsv(header: true)
            .map { row -> tuple(row.condition, file(row.fraction_file)) }
            .groupTuple(by: 0)
            .map { condition, fractions ->
                def meta = [condition: condition]
                tuple(meta, fractions)
            }

        by_count = rows.branch { meta, fs ->
            single: fs.size() == 1
            multi:  fs.size() > 1
        }

        MERGE_GPF(by_count.multi)

        single_ch = by_count.single.map { meta, fs -> tuple(meta, fs[0]) }
        combined  = MERGE_GPF.out.combined.mix(single_ch)

        SEARCH_GPF(
            combined.map { meta, c -> tuple(meta, c, fasta) },
            dlib_ch
        )

        pooled = SEARCH_GPF.out.products
            .map { meta, files -> files }
            .flatten()
            .collect()

        ENCYCLOPEDIA_LIBEXPORT(pooled, 'chromatogram', dlib_ch, fasta, params.gpf_search_args)
        library_ch = ENCYCLOPEDIA_LIBEXPORT.out.elib

    } else {
        library_ch = dlib_ch
    }

    emit:
    library                 = library_ch
    reference_list_derived  = reference_list_ch
}
