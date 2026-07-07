#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PREPARE_LIBRARY } from './subworkflows/local/prepare_library'
include { DIA_WORKFLOW    } from './subworkflows/local/dia_workflow'
include { PRM_WORKFLOW    } from './subworkflows/local/prm_workflow'
include { MSCONVERT       } from './modules/local/msconvert/msconvert'


def helpMessage() {
    log.info """
    Usage:
      nextflow run main.nf --input samplesheet.csv --fasta proteins.fasta [options]

    Samplesheet columns (CSV with header):
      sample_id   unique identifier
      file        path to raw MS file: .raw | .mzML | .dia
      mode        DIA | PRM
      mass_list   path to assay .txt/.csv (required for PRM, ignored for DIA)

    Library options (fastest path listed first):
      --elib            Pre-built chromatogram library (.elib) -- skips all library building
      --dlib            Pre-built spectral library (.dlib)
      --blib            Skyline .blib (converted to .dlib via -blibToLib)
      --use_koina       Predict .dlib from --fasta via Koina
      --library_sheet   CSV (columns: condition, fraction_file) for building an .elib
                        from GPF fractions. Multiple fractions per condition are merged
                        first, then searched. Requires one of --dlib/--blib/--use_koina
                        as the base .dlib.
    """.stripIndent()
}

workflow {
    if (!params.input || !params.fasta) { helpMessage(); exit 1 }

    samples = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [id: row.sample_id, mode: row.mode.toUpperCase()]
            def ml   = row.mass_list?.trim() ? file(row.mass_list) : file('NO_FILE')
            tuple(meta, file(row.file), ml)
        }

    fasta = file(params.fasta)

    PREPARE_LIBRARY(fasta)
    library = PREPARE_LIBRARY.out.library

    ext_split = samples.branch { meta, f, _ml ->
        raw:   f.getName().toLowerCase().endsWith('.raw')
        ready: true
    }

    MSCONVERT(ext_split.raw)

    normalized = ext_split.ready.mix(MSCONVERT.out.mzml)

    by_mode = normalized.branch { meta, _f, _ml ->
        dia: meta.mode == 'DIA'
        prm: meta.mode == 'PRM'
    }

    DIA_WORKFLOW(by_mode.dia, library, fasta)
    PRM_WORKFLOW(by_mode.prm, library, fasta)
}
