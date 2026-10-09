include { KOINA_PREDICT                } from '../../modules/local/encyclopedia/koina_predict'
include { BLIB_TO_DLIB                 } from '../../modules/local/encyclopedia/blib_to_dlib'
include { RESTRICT_BACKGROUND          } from '../../modules/local/library/restrict_background'
include { PREPARE_CHROMATOGRAM_LIBRARY } from './prepare_chromatogram_library'

workflow PREPARE_BACKGROUND {
    take:
    source
    background_library
    background_fasta

    main:
    def lib   = background_library ? file(background_library) : null
    def fasta = background_fasta   ? file(background_fasta)   : null
    def is_blib = lib && lib.name.toLowerCase().endsWith('.blib')
    def dlib_ch

    if (source == 'library') {
        if (is_blib) {
            BLIB_TO_DLIB(lib, fasta)
            dlib_ch = BLIB_TO_DLIB.out.dlib
        } else {
            dlib_ch = Channel.value(lib)
        }

    } else if (source == 'predicted') {
        KOINA_PREDICT(fasta)
        dlib_ch = KOINA_PREDICT.out.dlib
        if (params.background_peptides) {
            RESTRICT_BACKGROUND(dlib_ch, file(params.background_peptides),
                                params.background_peptides_top)
            dlib_ch = RESTRICT_BACKGROUND.out.dlib
        }

    } else if (source == 'measured') {
        def slib = params.background_search_library ? file(params.background_search_library) : null
        def search_lib
        if (slib && slib.name.toLowerCase().endsWith('.blib')) {
            BLIB_TO_DLIB(slib, fasta)
            search_lib = BLIB_TO_DLIB.out.dlib
        } else if (slib) {
            search_lib = Channel.value(slib)
        } else {
            KOINA_PREDICT(fasta)
            search_lib = KOINA_PREDICT.out.dlib
        }
        PREPARE_CHROMATOGRAM_LIBRARY(params.background_sheet, search_lib, fasta, 'background')
        dlib_ch = PREPARE_CHROMATOGRAM_LIBRARY.out.elib

    } else {
        error "PREPARE_BACKGROUND: unknown background_source '${source}'"
    }

    emit:
    dlib = dlib_ch
}
