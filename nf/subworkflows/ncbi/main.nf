#!/usr/bin/env nextflow
// main nextflow script for EGAPx execution
// route data to subworkflows

nextflow.enable.dsl=2
import groovy.yaml.YamlBuilder

include { rnaseq_short_plane } from './rnaseq_short/main'
include { rnaseq_long_plane } from './rnaseq_long/main'
include { target_proteins_plane } from './target_proteins/main'
include { gnomon_plane } from './gnomon/main'
include { orthology_plane } from './orthology/main'
include { annot_proc_plane } from './annot_proc/main'
include { setup_genome; setup_proteins } from './setup/main'
include { convert_annotations } from './default/convert_annotations/main'
include { convert_summary_files } from './default/convert_summary_files/main'
include { cmsearch_plane } from './cmsearch/main'
include { busco } from './busco/main'
include { winmask_plane } from './winmask/main'

params.verbose = false


workflow egapx {
    take:
        input_params     // input parameters
        // genome           // path to genome
        // proteins         // path to proteins, optional
        // proteins_trusted // path to trusted proteins, optional

        // // Alternative groups of parameters, one of them should be set
        // // reads_query - SRA query in the form accepted by NCBI
        // // reads_ids - list of SRA IDs
        // // reads, reads_metadata - path to reads accompanied by metadata
        // reads_query     // SRA query
        // reads_ids       // list of SRA IDs
        // reads           // path to reads
        // reads_metadata  // path to reads metadata 13 tab-delimited fields, 1-st - SRA ID, 3-rd paired or unpaired, everything else - not used, but must be present
        //                 // 4, 5, 13 - numbers, 5 - non zero number
        // organelles      // path to organelle list
        // // Alternative parameters, one of them should be set
        // // tax_id - NCBI tax id of the closest taxon to the genome
        // // hmm_params - HMM parameters
        // tax_id          // NCBI tax id of the closest taxon to the genome
        // symbol_format_class // string for how to create gene names
        // hmm_params      // HMM parameters
        // train_hmm       // Boolean, whether to train HMM
        // //
        // softmask        // softmask for GNOMON, optional
        // //
        // max_intron      // max intron length
        // genome_size_threshold // the threshold for calculating actual max intron length
        // ortho_files     // files reference genome sequence and annotation for find_orthology
        // reference_sets  // reference sets, for now only swissprot
        // prot_denylist   // path to protein denylist
        task_params     // task parameters for every task
    main:
        print "workflow.container: ${workflow.container}"

        def genome = input_params.get('genome', [])
        def proteins = input_params.get('proteins', [])
        def proteins_trusted = input_params.get('proteins_trusted', [])
        def short_reads_query = input_params.get('short_reads_query', [])
        def short_reads_ids = input_params.get('short_reads_ids', [])
        def short_reads = input_params.get('short_reads', [])
        def short_reads_metadata = input_params.get('short_reads_metadata', [])
        def long_reads_query = input_params.get('long_reads_query', [])
        def long_reads_ids = input_params.get('long_reads_ids', [])
        def long_reads = input_params.get('long_reads', [])
        def long_reads_metadata = input_params.get('long_reads_metadata', [])
        def organelles = input_params.get('organelles', []) ?: []
        def tax_id = input_params.get('taxid', [])
        def symbol_format_class = input_params.get('symbol_format_class', [])
        def name_cleanup_rules_file = input_params.get('name_cleanup_rules_file', [])
        def hmm_params = input_params.get('hmm', []) ?: []
        def train_hmm = input_params.get('train_hmm', [])
        def max_intron = input_params.get('max_intron', [])
        def genome_size_threshold = input_params.get('genome_size_threshold', [])
        def ortho_files = input_params.get('ortho', []) ?: []
        def rnaseq_alignments = input_params.get('rnaseq_alignments', []) ?: []
        def protein_alignments = input_params.get('protein_alignments', []) ?: []
        def reference_sets = input_params.get('reference_sets', []) ?: []
        def prot_denylist = input_params.get('prot_denylist', []) ?: []
        def export_bam = input_params.get('export_bam', []) ?: []
        def gnomon_filtering_scores_file = input_params.get('gnomon_filtering_scores_file', []) ?: []
        def locus_tag_prefix = input_params.get('locus_tag_prefix', []) ?: ''
        def busco_lineage = input_params.get('busco_lineage', '')
        def busco_lineage_download = input_params.get('busco_lineage_download', [])
        def annotation_provider = input_params.get('annotation_provider', []) ?: 'GenBank submitter'
        def annotation_name_prefix = input_params.get('annotation_name_prefix', []) ?: 'EGAPx Test Assembly'

        if (params.verbose) {
            // dump all params as yaml
            def builder = new YamlBuilder()
            builder(params)
            println(builder.toString())
        }
         
        setup_genome_params = task_params.get('setup', [:])
        setup_genome_params['max_intron'] = max_intron
        setup_genome_params['genome_size_threshold'] = genome_size_threshold
        setup_genome_params['annotation_provider'] = annotation_provider
        setup_genome_params['annotation_name_prefix'] = annotation_name_prefix
        (scaffolds, gencoll_asn, unpacked_genome, genome_asn, genome_asnb, eff_max_intron) = setup_genome(genome, organelles, setup_genome_params)

        // Protein alignments
        // def protein_alignments = []
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
        // def rnaseq_alignments = []
        def star_bam = []
        if (short_reads_query || short_reads_ids || short_reads) {
            rnaseq_short_plane(genome_asn, scaffolds, unpacked_genome, short_reads_query, short_reads_ids, short_reads, short_reads_metadata, organelles, tax_id, eff_max_intron, task_params) 
            rnaseq_short_alignments = rnaseq_short_plane.out.rnaseq_alignments
            sra_exons = rnaseq_short_plane.out.sra_exons
            sra_exons_slices = rnaseq_short_plane.out.sra_exons_slices
            if (export_bam){
                star_bam = rnaseq_short_plane.out.star_bam
            }
        }

        // Combine RNASeq short and protein alignments
        def rnsp_alignments
        if (proteins && (short_reads_query || short_reads_ids || short_reads)) [
            rnsp_alignments = rnaseq_short_alignments.combine(protein_alignments)
        ] else if (proteins) {
            rnsp_alignments = protein_alignments
        } else {
            rnsp_alignments = rnaseq_short_alignments
        }

        def alignments
        // RNASeq long alignments
        if (long_reads_query || long_reads_ids || long_reads) {
            rnaseq_long_plane(unpacked_genome, gencoll_asn, long_reads_query, long_reads_ids, long_reads, eff_max_intron, task_params)
            alignments = rnsp_alignments.combine(rnaseq_long_plane.out.alignments)
        } else {
            alignments = rnsp_alignments
        }

        // Winmask
        winmask_plane (genome_asnb, scaffolds, gencoll_asn, [], [], task_params)  
        win_softmask = winmask_plane.out.softmask

        // GNOMON

        def gnomon_models = []
        def effective_hmm
        gnomon_plane(genome_asn, scaffolds, gencoll_asn, proteins_asn, alignments, sra_exons, sra_exons_slices, proteins_trusted, tax_id, hmm_params, train_hmm, win_softmask, eff_max_intron, reference_sets, gnomon_filtering_scores_file, task_params) 
        gnomon_models = gnomon_plane.out.gnomon_models


        // TODO: enable and wire cmsearch_plane.out.cmsearch_annots into annot_builder
        // cmsearch_plane(unpacked_genome)

        // outputs 
       
        def accept_annot_file = []
        def gff_annotated_file = []
        def final_asn_out = []  
        def locus_out = []
        def stats_dir = []
        def annotated_genome_file = []
        def annotation_data_comment_file = []

        annot_proc_plane(annotation_name_prefix, gnomon_models, gencoll_asn, genome_asn, genome_asnb,
                         scaffolds, tax_id, symbol_format_class, name_cleanup_rules_file,
                         ortho_files, gnomon_plane.out.alignments, gnomon_plane.out.best_naming_hits, gnomon_plane.out.swiss_prot_asn, prot_denylist, task_params)

        locus_out = annot_proc_plane.out.locus
        locustypes = annot_proc_plane.out.locustypes
        final_asn_out = annot_proc_plane.out.final_asn_out
        accept_annot_file = annot_proc_plane.out.accept_annot_file
        gff_annotated_file = annot_proc_plane.out.gff_annotated_file
        stats_dir = annot_proc_plane.out.stats
        annotated_genome_file = annot_proc_plane.out.annotated_genome_asn
        annotation_data_comment_file = annot_proc_plane.out.annotation_data_comment

        convert_annotations(annot_proc_plane.out.to_convert, task_params.get('convert_annotations', [:])) 
        convert_summary_files(gnomon_plane.out.gnomon_report, gnomon_plane.out.gnomon_quality_report ,locustypes, locus_tag_prefix)

        // BUSCO
        def busco_out = []
        if (busco_lineage) {
            busco(annot_proc_plane.out.annot_proteins, busco_lineage, busco_lineage_download, task_params.get('busco', [:]))
            busco_out = busco.out.results
        } 

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
        out_star_bam = star_bam
        user_gnomon_report = convert_summary_files.out.gnomon_report
        user_gnomon_quality_report = convert_summary_files.out.gnomon_quality_report
        
        busco_results = busco_out
        gnomon_biotype_contam_rpt = annot_proc_plane.out.gnomon_biotype_contam_rpt 

        //converted_outs = converted_outs
}
