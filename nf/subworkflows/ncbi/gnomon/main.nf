#!/usr/bin/env nextflow
// gnomon plane workflow
// route data to tasks

nextflow.enable.dsl=2

params.import_prefix = "../../../../nf/subworkflows/ncbi/" // redirected during testing

include { chainer_wnode as chainer } from "./${params.import_prefix}gnomon/chainer_wnode/main"
include { gnomon_wnode } from "./${params.import_prefix}gnomon/gnomon_wnode/main"
include { gnomon_training_iterations } from "./${params.import_prefix}gnomon-training-iteration/gnomon_training_iterations/main"
include { gnomon_biotype} from "./${params.import_prefix}annot_proc/gnomon_biotype/main"
//include { locus_track } from './locus_track/main'
//include { locus_link } from './locus_link/main'


params.intermediate = false

workflow gnomon_plane {
    take:
        genome_asn
        scaffolds
        gencoll_asn
        proteins_asn
        alignments // list of all relevent input alignments

        proteins_trusted
        // Alternative parameters, one of them should be set
        // tax_id - NCBI tax id of the closest taxon to the genome
        // hmm_params - HMM parameters
        tax_id          // NCBI tax id of the closest taxon to the genome
        hmm_params      // HMM parameters
        train_hmm       // Boolean, whether to train HMM
        //
        softmask        // softmask for GNOMON, optional
        max_intron      // max intron length
        task_params     // task parameters for every task
    main:
        // GNOMON
        def effective_hmm
        if (!train_hmm) {
            effective_hmm = hmm_params
        } else {
            effective_hmm = gnomon_training_iterations(hmm_params, genome_asn, proteins_asn, alignments, /* evidence_denylist */ [], /* gap_fill_allowlist */ [],
                [proteins_trusted].flatten(), scaffolds, softmask, scaffolds,
                max_intron,
                task_params)
        }

        chainer(alignments, effective_hmm, /* evidence_denylist */ [], /* gap_fill_allowlist */ [], scaffolds, [proteins_trusted].flatten(), genome_asn, proteins_asn, task_params.get('chainer_wnode', [:]))

        def gn_models = []
        gnomon_wnode(scaffolds, chainer.out.chains, chainer.out.chains_slices, effective_hmm, [], softmask, genome_asn, proteins_asn, task_params.get('gnomon_wnode', [:]))
       
    emit:
        gnomon_models = gnomon_wnode.out.outputs
        // trained_hmm = effective_hmm
}



