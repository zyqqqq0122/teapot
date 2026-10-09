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
include { PREPARE_CHROMATOGRAM_LIBRARY            } from './prepare_chromatogram_library'
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
    def need_dlib   = params.library_sheet || !params.elib
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

                PREPARE_BACKGROUND(background_source, background_library, background_fasta)
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
        } else if (params.library_sheet) {
            error "--library_sheet needs a library to search its runs with: --dlib, or " +
                  "--use_koina to predict one from --fasta."
        } else {
            error "No spectral library available. Provide --elib, --dlib, --blib, " +
                  "--library_sheet or --use_koina."
        }
    }

    def library_ch
    if (params.library_sheet) {
        PREPARE_CHROMATOGRAM_LIBRARY(params.library_sheet, dlib_ch, fasta, 'chromatogram')
        library_ch = PREPARE_CHROMATOGRAM_LIBRARY.out.elib
    } else if (params.elib) {
        library_ch = Channel.value(file(params.elib))
    } else {
        library_ch = dlib_ch
    }

    emit:
    library                 = library_ch
    reference_list_derived  = reference_list_ch
}
