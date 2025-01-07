#!/usr/bin/env nextflow
// main nextflow script for EGAPx execution
// route data to subworkflows

nextflow.enable.dsl=2

include { rnaseq_short_plane } from './rnaseq_short/main'
include { target_proteins_plane } from './target_proteins/main'
include { gnomon_plane } from './gnomon/main'
include { orthology_plane } from './orthology/main'
include { annot_proc_plane } from './annot_proc/main'
include { setup_genome; setup_proteins } from './setup/main'
//include { annot_builder } from './annot_proc/annot_builder/main'
//include { final_asn_markup } from './annot_proc/final_asn/main'
//include { annotwriter } from './annot_proc/annotwriter/main'
include {convert_annotations } from './default/convert_annotations/main' 

params.intermediate = false
params.use_orthology = false
params.use_post_gnomon = false


workflow egapx {
    take:
        genome           // path to genome
        proteins         // path to proteins, optional
        proteins_trusted // path to trusted proteins, optional

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
        symbol_format_class // string for how to create gene names
        hmm_params      // HMM parameters
        train_hmm       // Boolean, whether to train HMM
        //
        softmask        // softmask for GNOMON, optional
        //
        max_intron      // max intron length
        genome_size_threshold // the threshold for calculating actual max intron length
        ortho_files     // files reference genome sequence and annotation for find_orthology
        reference_sets  // reference sets, for now only swissprot
        prot_denylist   // path to protein denylist
        task_params     // task parameters for every task
    main:
        print "workflow.container: ${workflow.container}"

        def setup_genome_params = task_params.get('setup', [:])
        setup_genome_params['max_intron'] = max_intron
        setup_genome_params['genome_size_threshold'] = genome_size_threshold
        (scaffolds, gencoll_asn, unpacked_genome, genome_asn, genome_asnb, eff_max_intron) = setup_genome(genome, organelles, setup_genome_params)

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
        gnomon_plane(genome_asn, scaffolds, gencoll_asn, proteins_asn, alignments, proteins_trusted, tax_id, hmm_params, train_hmm, softmask, eff_max_intron, task_params) 
        gnomon_models = gnomon_plane.out.gnomon_models


        // outputs 
       
        def accept_annot_file = []
        def gff_annotated_file = []
        def final_asn_out = []  
        def locus_out = []
        def stats_dir = []
        def annotated_genome_file = []
        def annotation_data_comment_file = []
        annot_proc_plane(gnomon_models, gencoll_asn, genome_asn, genome_asnb, scaffolds, tax_id, symbol_format_class, ortho_files, reference_sets, prot_denylist, task_params)
        locus_out = annot_proc_plane.out.locus
        final_asn_out = annot_proc_plane.out.final_asn_out
        accept_annot_file = annot_proc_plane.out.accept_annot_file
        gff_annotated_file = annot_proc_plane.out.gff_annotated_file
        stats_dir = annot_proc_plane.out.stats
        annotated_genome_file = annot_proc_plane.out.annotated_genome_asn
        annotation_data_comment_file = annot_proc_plane.out.annotation_data_comment

        convert_annotations(annot_proc_plane.out.to_convert, task_params.get('convert_annotations', [:])) 
        
    emit:
        out_files = gff_annotated_file
        out_ggff = convert_annotations.out.genomic_gff
        out_ggtf = convert_annotations.out.genomic_gtf
        out_gfa  = convert_annotations.out.genomic_fasta
        out_rna_fa = convert_annotations.out.transcripts_fasta
        out_cds_fa = convert_annotations.out.cds_fasta
        out_prot_fa = convert_annotations.out.proteins_fasta
        annot_builder_output = annot_proc_plane.out.accept_annot_file
        locus = locus_out
        final_asn_outputs = final_asn_out
        validated = annot_proc_plane.out.validated
        stats = stats_dir
        annotated_genome_asn = annotated_genome_file
        annotation_data_comment = annotation_data_comment_file
        gnomon_summaries = gnomon_plane.out.gnomon_summaries
        gnomon_quality_report = gnomon_plane.out.gnomon_quality_report
        gnomon_report = gnomon_plane.out.gnomon_report
        //converted_outs = converted_outs
}
