#!/usr/bin/env nextflow
// main nextflow script for EGAPx execution
// route data to subworkflows

nextflow.enable.dsl=2

include { rnaseq_short_plane } from './rnaseq_short/main'
include { target_proteins_plane } from './target_proteins/main'
include { gnomon_plane; post_gnomon_plane } from './gnomon/main'
include { orthology_plane } from './orthology/main'
include { setup_genome; setup_proteins } from './setup/main'
include { annot_builder } from './default/annot_builder/main'
include { annotwriter } from './default/annotwriter/main'


params.intermediate = false
params.use_orthology = false
params.use_post_gnomon = false


workflow egapx {
    take:
        genome          // path to genome
        proteins        // path to proteins, optional

        // Alternative groups of parameters, one of them should be set
        // reads_query - SRA query in the form accepted by NCBI
        // reads_ids - list of SRA IDs
        // reads, reads_metadata - path to reads accompanied by metadata
        reads_query     // SRA query
        reads_ids       // list of SRA IDs
        reads           // path to reads
        reads_metadata  // path to reads metadata 13 tab-delimited fields, 1-st - SRA ID, 3-rd paired or unpaired, everything else - not used, but must be present
                        // 4, 5, 13 - numbers, 5 - non zero number

        organelles      // path to organelle list
        // Alternative parameters, one of them should be set
        // tax_id - NCBI tax id of the closest taxon to the genome
        // hmm_params - HMM parameters
        tax_id          // NCBI tax id of the closest taxon to the genome
        hmm_params      // HMM parameters
        hmm_taxid       // NCBI tax id of the HMM
        //
        softmask        // softmask for GNOMON, optional
        //
        max_intron      // max intron length
        genome_size_threshold // the threshold for calculating actual max intron length
        task_params     // task parameters for every task
    main:
        print "workflow.container: ${workflow.container}"

        def setup_genome_params = task_params.get('setup', [:])
        setup_genome_params['max_intron'] = max_intron
        setup_genome_params['genome_size_threshold'] = genome_size_threshold
        def (scaffolds, gencoll_asn, unpacked_genome, genome_asn, genome_asnb, eff_max_intron) = setup_genome(genome, organelles, setup_genome_params)

        // Protein alignments
        def protein_alignments = []
        def unpacked_proteins
        def proteins_asn = []
        def proteins_asnb = []
        if (proteins) {
            // miniprot plane
            (unpacked_proteins, proteins_asn) = setup_proteins(proteins, task_params.get('setup', [:]))
            target_proteins_plane(unpacked_genome, genome_asn, gencoll_asn, unpacked_proteins, proteins_asn, eff_max_intron, task_params)
            protein_alignments = target_proteins_plane.out.protein_alignments
        }

        // RNASeq short alignments
        def rnaseq_alignments = []
        if (reads_query || reads_ids || reads) {
            rnaseq_short_plane(genome_asn, scaffolds, unpacked_genome, reads_query, reads_ids, reads, reads_metadata, organelles, tax_id, eff_max_intron, task_params) 
            rnaseq_alignments = rnaseq_short_plane.out.rnaseq_alignments
        }

        // Combine RNASeq and protein alignments
        def alignments
        if (proteins && (reads_query || reads_ids || reads)) [
            alignments = rnaseq_alignments.combine(protein_alignments)
        ] else if (proteins) {
            alignments = protein_alignments
        } else {
            alignments = rnaseq_alignments
        }

        // GNOMON

        def gnomon_models = []
        def effective_hmm
        gnomon_plane(genome_asn, scaffolds, gencoll_asn, proteins_asn, alignments, tax_id, hmm_params, hmm_taxid, softmask, eff_max_intron, task_params) 
        gnomon_models = gnomon_plane.out.gnomon_models


        // outputs 
        annot_builder(gencoll_asn, gnomon_models, genome_asn, task_params.get('annot_builder', [:]))
        def accept_annot_file = annot_builder.out.accept_ftable_annot 
        def annot_files = annot_builder.out.annot_files

        if (params.use_orthology) {
            // ORTHOLOGY
            orthology_plane(genome_asnb, gencoll_asn, gnomon_models, annot_files, task_params)
            def orthologs = orthology_plane.out.orthologs
            if (params.use_post_gnomon) {
                //POST GNOMON
                post_gnomon_plane(gnomon_models, gencoll_asn, orthologs, tax_id, task_params)
            }
        }

        annotwriter(accept_annot_file, [:])
        annotwriter.out.annoted_file     

    emit:
        out_files = annotwriter.out.annoted_file
        annot_builder_output = annot_builder.out.outputs
        // locus = post_gnomon_plane.out.locus
}
