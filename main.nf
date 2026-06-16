#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { OPENSWATH    } from './modules/openswath'
include { ENCYCLOPEDIA } from './modules/encyclopedia'
include { PYPROPHET    } from './modules/pyprophet'

workflow {

    // Validate required params
    if (!params.input)        error "Please provide --input (sample sheet CSV)"
    if (!params.spectral_lib) error "Please provide --spectral_lib"
    if (!params.fasta)        error "Please provide --fasta"

    // Parse sample sheet: columns — sample, mzml
    ch_samples = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> tuple(row.sample, file(row.mzml)) }

    ch_lib  = Channel.fromPath(params.spectral_lib)
    ch_fasta = Channel.fromPath(params.fasta)

    // Run tools
    OPENSWATH(ch_samples, ch_lib)
    ENCYCLOPEDIA(ch_samples, ch_lib, ch_fasta)
    PYPROPHET(OPENSWATH.out.osw)
}
