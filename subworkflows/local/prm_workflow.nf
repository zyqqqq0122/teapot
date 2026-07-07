include { CONTEXT_SEARCH } from '../../modules/local/encyclopedia/context_search'

workflow PRM_WORKFLOW {
    take:
    samples   // tuple(meta, prm_file, mass_list)
    elib
    fasta

    main:
    CONTEXT_SEARCH(samples, elib, fasta)

    emit:
    reference_targets    = CONTEXT_SEARCH.out.ref_targets
    reference_decoys     = CONTEXT_SEARCH.out.ref_decoys
    background_targets   = CONTEXT_SEARCH.out.bg_targets
    background_decoys    = CONTEXT_SEARCH.out.bg_decoys
    plots                = CONTEXT_SEARCH.out.plots
}
