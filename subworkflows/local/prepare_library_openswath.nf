include { BLIB_TO_DLIB                            } from '../../modules/local/encyclopedia/blib_to_dlib'
include { BLIB_TO_REFERENCE_LIST                  } from '../../modules/local/encyclopedia/blib_to_reference_list'
include { DLIB_TO_OPENSWATH_TSV                   } from '../../modules/local/encyclopedia/dlib_to_openswath_tsv'
include { COMPLETE_CHANNELS                        } from '../../modules/local/library/complete_channels'
include { EMIT_IRT_LIBRARY                        } from '../../modules/local/library/emit_irt_library'
include { DERIVE_IRT_FROM_REFERENCE_LIST          } from '../../modules/local/library/derive_irt_from_reference_list'
include { HARMONIZE_RT as HARMONIZE_RT_TARGET     } from '../../modules/local/library/harmonize_rt'
include { HARMONIZE_RT as HARMONIZE_RT_BACKGROUND } from '../../modules/local/library/harmonize_rt'
include { MERGE_LIBRARIES                           } from '../../modules/local/library/merge_libraries'
include { ASSERT_RT_COHERENT                      } from '../../modules/local/library/assert_rt_coherent'
include { ASSERT_BACKGROUND_ADEQUATE              } from '../../modules/local/library/assert_background_adequate'
include { OPENSWATH_ASSAY_GENERATOR               } from '../../modules/local/openms/openswath_assay_generator'
include { OPENSWATH_DECOY_GENERATOR               } from '../../modules/local/openms/openswath_decoy_generator'
include { PREPARE_BACKGROUND                      } from './prepare_background'
include { KOINA_PREDICT                           } from '../../modules/local/encyclopedia/koina_predict'
include { MERGE_GPF             as OPENSWATH_MERGE_GPF      } from '../../modules/local/encyclopedia/merge_gpf'
include { SEARCH_GPF            as OPENSWATH_SEARCH_GPF     } from '../../modules/local/encyclopedia/search_gpf'
include { ENCYCLOPEDIA_LIBEXPORT as OPENSWATH_GPF_LIBEXPORT } from '../../modules/local/encyclopedia/encyclopedia_libexport'

workflow PREPARE_LIBRARY_OPENSWATH {
    take:
    fasta
    heavy_label
    background_source
    background_library
    background_fasta
    background_min_targets
    reference_list_for_anchors
    no_file

    main:
    def unimod_ch     = Channel.value(file("${projectDir}/assets/unimod_delta_masses.tsv"))
    def no_file_ch    = Channel.value(no_file)
    def run_channel   = (heavy_label as List)?.size() > 0
    def has_bg        = background_source && background_source != 'none'

    def library_ch
    def reference_list_ch = no_file_ch
    def irt_library_ch    = no_file_ch
    def full_window_ch    = no_file_ch
    def enc_lib_ch        = no_file_ch

    if (params.openswath_pqp) {
        library_ch = Channel.value(file(params.openswath_pqp))
        if (params.irt_library) {
            irt_library_ch = Channel.value(file(params.irt_library))
        }

    } else {
        def tsv_ch
        if (params.openswath_tsv || params.openswath_traml) {
            def osw_native = params.openswath_tsv ?: params.openswath_traml
            if (run_channel) {
                log.warn "PREPARE_LIBRARY_OPENSWATH: heavy_label is set but --openswath_tsv/--openswath_traml " +
                         "passes through as-is. COMPLETE_CHANNELS only operates on " +
                         "blib/dlib libraries; ensure your TSV contains both channels."
            }
            if (has_bg) {
                log.warn "PREPARE_LIBRARY_OPENSWATH: --openswath_tsv/--openswath_traml passes through; " +
                         "MERGE_LIBRARIES and HARMONIZE_RT are dlib/blib-only. " +
                         "Ensure your TSV already contains the background union " +
                         "in a coherent RT space."
            }
            tsv_ch = Channel.value(file(osw_native))
            if (params.irt_library) {
                irt_library_ch = Channel.value(file(params.irt_library))
            }

        } else if (params.openswath_blib || params.blib) {
            def blib_ch = Channel.value(file(params.openswath_blib ?: params.blib))

            def anchors_ch
            if (params.irt_library) {
                irt_library_ch = Channel.value(file(params.irt_library))
                anchors_ch     = params.irt_anchors ?
                                    Channel.value(file(params.irt_anchors)) :
                                    no_file_ch
                if (anchors_ch.equals(no_file_ch)) {
                    log.warn "PREPARE_LIBRARY_OPENSWATH: --irt_library supplied but " +
                             "--irt_anchors is not. HARMONIZE_RT will not have " +
                             "anchors to fit against. Set --irt_anchors or let " +
                             "the pipeline derive both from the vendor blib."
                }
            } else {
                EMIT_IRT_LIBRARY(blib_ch, heavy_label)
                anchors_ch     = EMIT_IRT_LIBRARY.out.anchors
                irt_library_ch = EMIT_IRT_LIBRARY.out.osw_tsv
            }

            def source_blib_ch
            if (run_channel) {
                COMPLETE_CHANNELS(blib_ch, heavy_label)
                source_blib_ch = COMPLETE_CHANNELS.out.library
            } else {
                source_blib_ch = blib_ch
            }
            BLIB_TO_DLIB(source_blib_ch, fasta)

            BLIB_TO_REFERENCE_LIST(source_blib_ch)
            reference_list_ch = BLIB_TO_REFERENCE_LIST.out.reference_list

            HARMONIZE_RT_TARGET(BLIB_TO_DLIB.out.dlib, anchors_ch, 'target',
                                false, false)

            def dlib_ch
            if (has_bg) {
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

                full_window_ch = HARMONIZE_RT_TARGET.out.full_window
                    .mix(HARMONIZE_RT_BACKGROUND.out.full_window)
                    .first()
                    .ifEmpty(no_file)
            } else {
                if (run_channel) {
                    log.warn "PREPARE_LIBRARY_OPENSWATH: --openswath_blib + heavy_label " +
                             "without a background library will produce " +
                             "reference_list == library and zero background " +
                             "targets in SPLIT_OPENSWATH_FEATURES. Set background_source " +
                             "to 'koina' | 'dlib' | 'blib' | 'elib' | 'gpf'."
                }
                dlib_ch = HARMONIZE_RT_TARGET.out.library
                full_window_ch = HARMONIZE_RT_TARGET.out.full_window
                    .ifEmpty(no_file)
            }

            enc_lib_ch = dlib_ch
            DLIB_TO_OPENSWATH_TSV(dlib_ch, unimod_ch)
            tsv_ch = DLIB_TO_OPENSWATH_TSV.out.tsv

        } else if (params.use_koina || params.dlib || params.elib ||
                   (params.library_sheet && background_source != 'gpf')) {
            def has_user_irt      = params.irt_library && params.irt_anchors
            def has_user_only_lib = params.irt_library && !params.irt_anchors
            def has_ref_anchors   = reference_list_for_anchors &&
                                    reference_list_for_anchors.name != 'NO_FILE'
            def use_full_window   = params.rt_harmonize_fallback == 'full_window'

            if (has_user_only_lib) {
                error """
                    PREPARE_LIBRARY_OPENSWATH: --use_koina + --irt_library requires
                    --irt_anchors alongside it. The anchors TSV (columns:
                    PeptideModSeq, Irt, [Standard]) defines the reference scale
                    HARMONIZE_RT regresses the Koina library onto. Without it
                    the pipeline would fall back to passing Koina-iRT through
                    unregressed. If you reused a vendor blib's anchor set
                    (which often has its own custom iRT scale), that would
                    silently miscalibrate. Supply --irt_anchors, or drop
                    --irt_library and let the pipeline derive anchors from
                    the samplesheet's reference_list (Case-1 scheduled-PRM
                    default), or use the full-window fallback.
                    """.stripIndent().trim()
            }
            if (!has_user_irt && !has_ref_anchors && !use_full_window) {
                def route = params.use_koina ? '--use_koina'
                          : params.dlib      ? '--dlib'
                          : params.elib      ? '--elib' : '--library_sheet'
                error """
                    PREPARE_LIBRARY_OPENSWATH: ${route} cannot proceed without an
                    anchor source. A predicted or pre-built library's RT column is in
                    Prosit-Biognosys iRT units (~-50..150), not gradient
                    minutes; running OpenSwathWorkflow without an iRT
                    normalisation library would extract every peak group at
                    a nonsense retention time and pick essentially nothing.

                    Choose one:
                      * supply --irt_library <TraML|TSV|PQP> AND
                        --irt_anchors <stripped_seq<TAB>iRT TSV>;
                      * put a reference_list with an RT column on your PRM
                        samples in the samplesheet. The pipeline derives
                        anchors automatically (Case-1 default; works whenever
                        the scheduled PRM method has RT windows);
                      * set --rt_harmonize_fallback full_window to accept
                        -rt_extraction_window -1 (correct but expensive per
                        sample; loud warning on large libraries).
                    """.stripIndent().trim()
            }

            def src_lib_ch
            if (params.library_sheet && background_source != 'gpf') {
                def seed_ch
                if (params.dlib) {
                    seed_ch = Channel.value(file(params.dlib))
                } else if (params.use_koina) {
                    KOINA_PREDICT(fasta)
                    seed_ch = KOINA_PREDICT.out.dlib
                } else {
                    error "PREPARE_LIBRARY_OPENSWATH: --library_sheet needs a seed " +
                          "spectral library. Supply --dlib or --use_koina."
                }
                gpf_rows = Channel
                    .fromPath(params.library_sheet)
                    .splitCsv(header: true)
                    .map { row -> tuple(row.condition, file(row.fraction_file)) }
                    .groupTuple(by: 0)
                    .map { condition, fractions -> tuple([condition: condition], fractions) }
                gpf_split = gpf_rows.branch { _meta, fs ->
                    single: fs.size() == 1
                    multi:  fs.size() > 1
                }
                OPENSWATH_MERGE_GPF(gpf_split.multi)
                gpf_combined = OPENSWATH_MERGE_GPF.out.combined
                    .mix(gpf_split.single.map { meta, fs -> tuple(meta, fs[0]) })
                OPENSWATH_SEARCH_GPF(gpf_combined.map { meta, c -> tuple(meta, c, fasta) }, seed_ch)
                gpf_pooled = OPENSWATH_SEARCH_GPF.out.products
                    .map { _meta, files -> files }
                    .flatten()
                    .collect()
                OPENSWATH_GPF_LIBEXPORT(gpf_pooled, 'chromatogram_osw', seed_ch, fasta, params.gpf_search_args)
                src_lib_ch = OPENSWATH_GPF_LIBEXPORT.out.elib

            } else if (params.elib) {
                src_lib_ch = Channel.value(file(params.elib))
            } else if (params.dlib) {
                src_lib_ch = Channel.value(file(params.dlib))
            } else {
                KOINA_PREDICT(fasta)
                src_lib_ch = KOINA_PREDICT.out.dlib
            }

            def source_dlib_ch
            if (run_channel) {
                COMPLETE_CHANNELS(src_lib_ch, heavy_label)
                source_dlib_ch = COMPLETE_CHANNELS.out.library
            } else {
                source_dlib_ch = src_lib_ch
            }

            def anchors2_ch = no_file_ch
            if (has_user_irt) {
                anchors2_ch    = Channel.value(file(params.irt_anchors))
                irt_library_ch = Channel.value(file(params.irt_library))
                HARMONIZE_RT_TARGET(source_dlib_ch, anchors2_ch,
                                    'target', false, false)
            } else if (has_ref_anchors) {
                DERIVE_IRT_FROM_REFERENCE_LIST(reference_list_for_anchors, src_lib_ch)
                anchors2_ch    = DERIVE_IRT_FROM_REFERENCE_LIST.out.anchors
                irt_library_ch = DERIVE_IRT_FROM_REFERENCE_LIST.out.osw_tsv
                HARMONIZE_RT_TARGET(source_dlib_ch, anchors2_ch,
                                    'target', false, false)
            } else {
                HARMONIZE_RT_TARGET(source_dlib_ch, no_file_ch, 'target',
                                    true, true)
            }

            def dlib2_ch
            if (has_bg) {
                PREPARE_BACKGROUND(fasta, background_source,
                                   background_library, background_fasta)
                HARMONIZE_RT_BACKGROUND(PREPARE_BACKGROUND.out.dlib,
                                        anchors2_ch, 'background',
                                        false, false)
                MERGE_LIBRARIES(HARMONIZE_RT_TARGET.out.library,
                              HARMONIZE_RT_BACKGROUND.out.library)
                dlib2_ch = MERGE_LIBRARIES.out.dlib
                ASSERT_RT_COHERENT(dlib2_ch)
                if (has_ref_anchors) {
                    ASSERT_BACKGROUND_ADEQUATE(dlib2_ch,
                                               Channel.value(reference_list_for_anchors),
                                               background_min_targets)
                }
            } else {
                dlib2_ch = HARMONIZE_RT_TARGET.out.library
            }

            enc_lib_ch = dlib2_ch
            DLIB_TO_OPENSWATH_TSV(dlib2_ch, unimod_ch)
            tsv_ch = DLIB_TO_OPENSWATH_TSV.out.tsv

            full_window_ch = HARMONIZE_RT_TARGET.out.full_window
                .ifEmpty(no_file)

        } else {
            error "No OpenSWATH library source. Provide an OSW-native library " +
                  "(--openswath_pqp, --openswath_tsv, --openswath_traml) or any engine-agnostic " +
                  "spectral-library source (--openswath_blib/--blib, --dlib, " +
                  "--elib, --library_sheet, --use_koina)."
        }

        OPENSWATH_ASSAY_GENERATOR(tsv_ch)
        OPENSWATH_DECOY_GENERATOR(OPENSWATH_ASSAY_GENERATOR.out.pqp)
        library_ch = OPENSWATH_DECOY_GENERATOR.out.pqp
    }

    emit:
    library                 = library_ch
    reference_list_derived  = reference_list_ch
    irt_library             = irt_library_ch
    full_window_flag        = full_window_ch
    encyclopedia_library    = enc_lib_ch
}
