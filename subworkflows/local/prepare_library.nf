include { KOINA_PREDICT } from '../../modules/local/encyclopedia/koina_predict'
include { BLIB_TO_DLIB  } from '../../modules/local/encyclopedia/blib_to_dlib'
include { MERGE_GPF     } from '../../modules/local/encyclopedia/merge_gpf'
include { SEARCH_GPF    } from '../../modules/local/encyclopedia/search_gpf'
include { BUILD_CHRLIB  } from '../../modules/local/encyclopedia/build_chrlib'

workflow PREPARE_LIBRARY {
    take:
    fasta

    main:
    def need_dlib = (params.library_sheet as boolean) || !params.elib

    def dlib_ch = null
    if (need_dlib) {
        if (params.dlib) {
            dlib_ch = Channel.value(file(params.dlib))
        } else if (params.blib) {
            BLIB_TO_DLIB(file(params.blib), fasta)
            dlib_ch = BLIB_TO_DLIB.out.dlib
        } else if (params.use_koina) {
            KOINA_PREDICT(fasta)
            dlib_ch = KOINA_PREDICT.out.dlib
        } else {
            error "No spectral library available. Provide --dlib, --blib, --use_koina, " +
                  "or supply --elib to skip .dlib entirely when not building from GPF."
        }
    }

    def library_ch
    if (params.elib && !params.library_sheet) {
        library_ch = Channel.value(file(params.elib))

    } else if (params.library_sheet) {
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

        BUILD_CHRLIB(pooled, dlib_ch, fasta)
        library_ch = BUILD_CHRLIB.out.elib

    } else {
        library_ch = dlib_ch
    }

    emit:
    library = library_ch
}
