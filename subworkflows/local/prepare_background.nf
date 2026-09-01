include { KOINA_PREDICT          } from '../../modules/local/encyclopedia/koina_predict'
include { BLIB_TO_DLIB           } from '../../modules/local/encyclopedia/blib_to_dlib'
include { MERGE_GPF              } from '../../modules/local/encyclopedia/merge_gpf'
include { SEARCH_GPF             } from '../../modules/local/encyclopedia/search_gpf'
include { ENCYCLOPEDIA_LIBEXPORT } from '../../modules/local/encyclopedia/encyclopedia_libexport'

workflow PREPARE_BACKGROUND {
    take:
    fasta
    source
    background_library_path
    background_fasta_path

    main:
    def dlib_ch = null

    if (source == 'none') {
        dlib_ch = null

    } else if (source == 'koina') {
        if (!background_fasta_path) {
            error "background_source='koina' requires --background_fasta. " +
                  "Do NOT reuse --fasta if it is a small vendor construct FASTA. " +
                  "Predicting a 'background' from the vendor constructs themselves " +
                  "yields a background composed entirely of vendor targets, which " +
                  "collides with the vendor library and gives zero background " +
                  "targets. Supply a proteome-wide FASTA."
        }
        def bg_fasta_ch = Channel.value(file(background_fasta_path))
        KOINA_PREDICT(bg_fasta_ch)
        dlib_ch = KOINA_PREDICT.out.dlib

    } else if (source == 'dlib') {
        if (!background_library_path) {
            error "background_source='dlib' requires --background_library."
        }
        dlib_ch = Channel.value(file(background_library_path))

    } else if (source == 'blib') {
        if (!background_library_path) {
            error "background_source='blib' requires --background_library."
        }
        BLIB_TO_DLIB(file(background_library_path), fasta)
        dlib_ch = BLIB_TO_DLIB.out.dlib

    } else if (source == 'elib') {
        if (!background_library_path) {
            error "background_source='elib' requires --background_library."
        }
        dlib_ch = Channel.value(file(background_library_path))

    } else if (source == 'gpf') {
        if (!params.library_sheet) {
            error "background_source='gpf' requires --library_sheet " +
                  "listing the GPF fractions to search."
        }

        def seed_dlib_ch
        if (params.dlib) {
            seed_dlib_ch = Channel.value(file(params.dlib))
            if (params.blib || params.openswath_blib) {
                log.warn "background_source='gpf' with --dlib as seed and " +
                         "vendor --blib/--openswath_blib as target: OK. " +
                         "Just verify --dlib is a proteome-wide seed and " +
                         "not derived from the vendor blib itself. An " +
                         "unspiked GPF run would filter the vendor heavies " +
                         "back out."
            }
        } else if (params.use_koina) {
            if (!background_fasta_path) {
                error "background_source='gpf' + --use_koina requires " +
                      "--background_fasta (proteome-wide). Reusing --fasta " +
                      "is unsafe if it is a small vendor FASTA."
            }
            def seed_fasta_ch = Channel.value(file(background_fasta_path))
            KOINA_PREDICT(seed_fasta_ch)
            seed_dlib_ch = KOINA_PREDICT.out.dlib
        } else {
            error "background_source='gpf' needs a seed library: either " +
                  "--dlib or --use_koina + --background_fasta."
        }

        def rows = Channel
            .fromPath(params.library_sheet)
            .splitCsv(header: true)
            .map { row -> tuple(row.condition, file(row.fraction_file)) }
            .groupTuple(by: 0)
            .map { condition, fractions ->
                tuple([condition: condition], fractions)
            }

        def by_count = rows.branch { _meta, fs ->
            single: fs.size() == 1
            multi:  fs.size() > 1
        }

        MERGE_GPF(by_count.multi)
        def single_ch = by_count.single.map { meta, fs -> tuple(meta, fs[0]) }
        def combined  = MERGE_GPF.out.combined.mix(single_ch)

        SEARCH_GPF(
            combined.map { meta, c -> tuple(meta, c, fasta) },
            seed_dlib_ch
        )

        def pooled = SEARCH_GPF.out.products
            .map { _meta, files -> files }
            .flatten()
            .collect()

        ENCYCLOPEDIA_LIBEXPORT(pooled, 'background_gpf',
                               seed_dlib_ch, fasta, params.gpf_search_args)
        dlib_ch = ENCYCLOPEDIA_LIBEXPORT.out.elib

    } else {
        error "Unknown background_source: '${source}'. " +
              "Valid values: none, koina, dlib, blib, elib, gpf."
    }

    emit:
    dlib = dlib_ch
}
