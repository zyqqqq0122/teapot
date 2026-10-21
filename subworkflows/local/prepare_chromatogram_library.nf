include { MSCONVERT              } from '../../modules/local/msconvert/msconvert'
include { MERGE_GPF              } from '../../modules/local/encyclopedia/merge_gpf'
include { SEARCH_GPF             } from '../../modules/local/encyclopedia/search_gpf'
include { ENCYCLOPEDIA_LIBEXPORT } from '../../modules/local/encyclopedia/encyclopedia_libexport'
include { CALIBRATE_FRAGMENT_MZ  } from '../../modules/local/runs/calibrate_fragment_mz'
include { DERIVE_SEARCH_TOLERANCE } from '../../modules/local/runs/derive_search_tolerance'

workflow PREPARE_CHROMATOGRAM_LIBRARY {
    take:
    sheet
    library
    fasta
    name
    calibrate

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

    def run_files = runs.ready.mix(converted)
    def warn_ch = Channel.empty()
    def tol_ch = Channel.value([value: null, runs: [:]])
    if (calibrate && Args.mode(params.fragment_mz_calibration) == 'auto') {
        def no_file = file("${projectDir}/assets/NO_FILE")
        CALIBRATE_FRAGMENT_MZ(run_files.map { condition, f ->
                                  tuple([id: f.baseName, condition: condition], f, no_file, no_file,
                                        Args.fragmentTolerance(params.sheet_search_args)) },
                              library, true)
        CALIBRATE_FRAGMENT_MZ.out.row
            .collectFile(name: "fragment_mz_calibration.${name}.tsv", storeDir: "${params.outdir}/qc",
                         keepHeader: true, skip: 1, sort: true)
        warn_ch = CALIBRATE_FRAGMENT_MZ.out.warnings.flatMap { Warnings.inFile(it, "calibration (${name})") }
        run_files = CALIBRATE_FRAGMENT_MZ.out.run.map { meta, f, _rl -> tuple(meta.condition, f) }
        DERIVE_SEARCH_TOLERANCE(CALIBRATE_FRAGMENT_MZ.out.tolerance.map { _m, f -> f }.collect().ifEmpty([no_file]),
                                no_file, Args.fragmentTolerance(params.sheet_search_args), name)
        tol_ch = DERIVE_SEARCH_TOLERANCE.out.tolerance.combine(DERIVE_SEARCH_TOLERANCE.out.runs)
            .map { t, r -> Args.readTolerance(t, r) }
            .first()
        tol_ch.subscribe { t -> log.info Warnings.fragmentTolerance(t, name, params.sheet_search_args) }
    }
    def args_of = { Map t, List fs -> fs.every { t.runs.containsKey(it.baseName) } ? Args.withBatchTolerance(params.sheet_search_args, t)
                                                                                : Args.flat(params.sheet_search_args) }

    by_count = run_files
        .groupTuple(by: 0)
        .combine(tol_ch)
        .map { condition, fractions, t -> tuple([condition: condition], fractions, args_of(t, fractions)) }
        .branch { _meta, fs, _args ->
            single: fs.size() == 1
            multi:  fs.size() > 1
        }

    MERGE_GPF(by_count.multi.map { meta, fs, _args -> tuple(meta, fs) })
    combined = MERGE_GPF.out.combined.join(by_count.multi.map { meta, _fs, args -> tuple(meta, args) })
        .mix(by_count.single.map { meta, fs, args -> tuple(meta, fs[0], args) })

    SEARCH_GPF(combined.map { meta, c, args -> tuple(meta, c, fasta, args) }, library)

    pooled = SEARCH_GPF.out.products
        .map { _meta, files -> files }
        .flatten()
        .collect()

    ENCYCLOPEDIA_LIBEXPORT(pooled, name, library, fasta, tol_ch.map { t -> Args.withBatchTolerance(params.sheet_search_args, t) }, false)

    emit:
    elib     = ENCYCLOPEDIA_LIBEXPORT.out.elib
    warnings = warn_ch
}
