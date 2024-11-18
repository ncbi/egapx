#!/usr/bin/env nextflow
// gnomon-only nextflow script for EGAPx execution
// route data to subworkflows

nextflow.enable.dsl=2

include { setup_genome; setup_proteins } from './setup/main'
include { get_hmm_params; run_get_hmm } from './default/get_hmm_params/main'
include { chainer_wnode as chainer } from './gnomon/chainer_wnode/main'
include { gnomon_wnode } from './gnomon/gnomon_wnode/main'
include { prot_gnomon_prepare } from './annot_proc/prot_gnomon_prepare/main'
include { annot_builder } from './annot_proc/annot_builder/main'
include { annotwriter } from './annot_proc/annotwriter/main'
include { run_align_sort} from './default/align_sort_sa/main'

params.intermediate = false

workflow only_gnomon {
    take:
        genome          // path to genome
        proteins        // path to proteins, optional
        // Alternative groups of parameters, one of them should be set
        rnaseq_alignments // path to rnaseq_collapse'ed alignments
        protein_alignments // path to miniprot, filtered, sorted, alignments
        
        organelles      // path to organelle list
        // Alternative parameters, one of them should be set
        // tax_id - NCBI tax id of the closest taxon to the genome
        // hmm_params - HMM parameters
        tax_id          // NCBI tax id of the closest taxon to the genome
        hmm_params      // HMM parameters
        train_hmm       // Boolean, whether to train HMM
        //
        softmask        // softmask for GNOMON, optional
        task_params     // task parameters for every task
    main:
        
        def (scaffolds, gencoll_asn, unpacked_genome, genome_asn) = setup_genome(genome, organelles, task_params.get('setup', [:]))

        // Protein alignments
        def unpacked_proteins
        def proteins_asn = []
        if (proteins) {
            // miniprot plane
            (unpacked_proteins, proteins_asn) = setup_proteins(proteins, task_params.get('setup', [:]))
        }

        // Combine RNASeq and protein alignments

        def alignments
        if (protein_alignments && rnaseq_alignments) {
            print(rnaseq_alignments.getClass())
            print(rnaseq_alignments)
            print(protein_alignments.getClass())
            print(protein_alignments)
            alignments = Channel.of(rnaseq_alignments).combine(Channel.of(protein_alignments))
        } else if (protein_alignments) {
            alignments = protein_alignments
        } else if (rnaseq_alignments) {
            alignments = rnaseq_alignments
        } else {
            print("error") 
        }

        // GNOMON

        def effective_hmm
        if (hmm_params) {
            effective_hmm = hmm_params
        } else {
            tmp_hmm = run_get_hmm(tax_id)
            b = tmp_hmm | splitText( { it.split('\n') } ) | flatten 
            c = b | last
            effective_hmm = c
        }

        chainer(alignments, effective_hmm, /* evidence_denylist */ [], /* gap_fill_allowlist */ [], scaffolds, /* trusted_genes */ [], genome_asn, proteins_asn, task_params.get('chainer', [:]))

        gnomon_wnode(scaffolds, chainer.out.chains, chainer.out.chains_slices, effective_hmm, [], softmask, genome_asn, proteins_asn, task_params.get('gnomon', [:]))
        def models = gnomon_wnode.out.outputs

        // prot_gnomon_prepare(models, task_params.get('prot_gnomon_prepare', [:]))
        
        // actual gnomon end but whatever


        annot_builder(gencoll_asn, models, genome_asn, task_params.get('annot_builder', [:]))
        def accept_asn = annot_builder.out.accept_asn
        
        annotwriter(accept_asn, [:])
        annotwriter.out.annoted_file
    emit:
        out_files = annotwriter.out.annoted_file
        evidence = annot_builder.out.outputs
}
