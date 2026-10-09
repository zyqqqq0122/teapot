include { MSCONVERT              } from '../../modules/local/msconvert/msconvert'
include { MERGE_GPF              } from '../../modules/local/encyclopedia/merge_gpf'
include { SEARCH_GPF             } from '../../modules/local/encyclopedia/search_gpf'
include { ENCYCLOPEDIA_LIBEXPORT } from '../../modules/local/encyclopedia/encyclopedia_libexport'

workflow PREPARE_CHROMATOGRAM_LIBRARY {
    take:
    sheet
    library
    fasta
    name

    main:
    def vendor_ext  = ['.raw', '.wiff', '.wiff2', '.d', '.lcd']
    def is_vendor   = { f -> vendor_ext.any { f.getName().toLowerCase().endsWith(it) } }
    def sidecars_of = { f ->
        ['.scan', '2.scan', '.timeseries.data'].collect { file("${f}${it}") }
                                               .findAll { it.exists() }
    }

    runs = Channel
        .fromPath(sheet)
        .splitCsv(header: true)
        .map { row -> tuple(row.condition, file(row.fraction_file)) }
        .branch { _condition, f ->
            vendor: is_vendor(f)
            ready:  true
        }

    MSCONVERT(runs.vendor.map { condition, f ->
        tuple([id: f.baseName, condition: condition], f, sidecars_of(f),
              file("${projectDir}/assets/NO_FILE"))
    })
    converted = MSCONVERT.out.mzml.flatMap { meta, mzml, _rl ->
        (mzml instanceof List ? mzml : [mzml]).collect { tuple(meta.condition, it) }
    }

    by_count = runs.ready
        .mix(converted)
        .groupTuple(by: 0)
        .map { condition, fractions -> tuple([condition: condition], fractions) }
        .branch { _meta, fs ->
            single: fs.size() == 1
            multi:  fs.size() > 1
        }

    MERGE_GPF(by_count.multi)
    combined = MERGE_GPF.out.combined.mix(by_count.single.map { meta, fs -> tuple(meta, fs[0]) })

    SEARCH_GPF(combined.map { meta, c -> tuple(meta, c, fasta) }, library)

    pooled = SEARCH_GPF.out.products
        .map { _meta, files -> files }
        .flatten()
        .collect()

    ENCYCLOPEDIA_LIBEXPORT(pooled, name, library, fasta, params.sheet_search_args, false)

    emit:
    elib = ENCYCLOPEDIA_LIBEXPORT.out.elib
}
