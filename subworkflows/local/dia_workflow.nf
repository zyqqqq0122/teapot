include { SEARCH_DIA } from '../../modules/local/encyclopedia/search_dia'

workflow DIA_WORKFLOW {
    take:
    samples   // tuple(meta, dia_or_mzml, _mass_list)
    library
    fasta

    main:
    dia_in = samples.map { meta, dia, _ml -> tuple(meta, dia) }
    SEARCH_DIA(dia_in, library, fasta)

    emit:
    features = SEARCH_DIA.out.features   // raw scored features
    targets  = SEARCH_DIA.out.targets    // FDR-filtered target peptides
    decoys   = SEARCH_DIA.out.decoys     // FDR-filtered decoys (for QC)
    elib     = SEARCH_DIA.out.elib       // per-file elib output
}
