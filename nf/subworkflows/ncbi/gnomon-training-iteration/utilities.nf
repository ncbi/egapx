#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// nextflow.preview.recursion=true

include { chainer_wnode as chainer } from '../gnomon/chainer_wnode/main'
include { gnomon_wnode } from '../gnomon/gnomon_wnode/main'
include { gnomon_training } from '../gnomon/gnomon_training/main'


workflow gnomon_training_iteration {
    take:
        initial_hmm_params
        genome_asn
        proteins_asn
        chainer_alignments
        chainer_evidence_denylist
        chainer_gap_fill_allowlist
        chainer_trusted_genes
        chainer_scaffolds
        gnomon_softmask
        gnomon_scaffolds
        max_intron
        parameters
    main:

        chainer(chainer_alignments, initial_hmm_params, chainer_evidence_denylist, chainer_gap_fill_allowlist, chainer_scaffolds, chainer_trusted_genes, genome_asn, proteins_asn, parameters.get('chainer_wnode', [:]))
        gnomon_wnode(gnomon_scaffolds, chainer.out.chains, chainer.out.chains_slices, initial_hmm_params, gnomon_softmask, [], genome_asn, proteins_asn,  parameters.get('gnomon_wnode', [:]))
        gnomon_training(genome_asn, gnomon_wnode.out.gn_models, max_intron, parameters.get('gnomon_training', [:]))

    emit:
        hmm_params_file = gnomon_training.out.hmm_params_file
        genome_asn = genome_asn
        proteins_asn = proteins_asn
        chainer_alignments = chainer_alignments
        chainer_evidence_denylist = chainer_evidence_denylist
        chainer_gap_fill_allowlist = chainer_gap_fill_allowlist
        chainer_trusted_genes = chainer_trusted_genes
        chainer_scaffolds = chainer_scaffolds
        gnomon_softmask = gnomon_softmask
        gnomon_scaffolds = gnomon_scaffolds
        max_intron = max_intron
        parameters = parameters
}
