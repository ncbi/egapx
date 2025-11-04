#!/usr/bin/env nextflow
// main nextflow script for EGAPx execution
// route data to subworkflows

nextflow.enable.dsl=2

params.import_prefix = "../../../../nf/subworkflows/ncbi/" // redirected during testing

include { miniprot } from "./${params.import_prefix}target_proteins/miniprot/main"
include { paf2asn } from "./${params.import_prefix}target_proteins/paf2asn/main"
include { tblastn_align } from "./${params.import_prefix}target_proteins/tblastn_align/main"
include { prosplign_prepare } from "./${params.import_prefix}target_proteins/prosplign_prepare/main"
include { prosplign_wnode } from "./${params.import_prefix}target_proteins/prosplign_wnode/main"
include { best_aligned_prot } from "./${params.import_prefix}target_proteins/best_aligned_prot/main"
include { align_filter_sa } from "./${params.import_prefix}target_proteins/align_filter_sa/main"
include { align_sort_sa} from "./${params.import_prefix}target_proteins/../default/align_sort_sa/main"


params.intermediate = false


workflow target_proteins_plane {
    take:
        unpacked_genome_fasta
        genome_asn
        genome_blastdb
        gencoll_asn
        unpacked_proteins_fasta
        proteins_asn
        aligner_name
        max_intron
        task_params     // task parameters for every task
    main:
        // Protein alignments
       
        alignments_to_use = []
        if (aligner_name == "miniprot" ) {
            miniprot(unpacked_genome_fasta, unpacked_proteins_fasta, max_intron, task_params.get('miniprot', [:]))
            def miniprot_file = miniprot.out.miniprot_file
            paf2asn(genome_asn, proteins_asn, miniprot_file, task_params.get('paf2asn', [:]))
            alignments_to_use = paf2asn.out.asn_file.collect()
        } 
        else if (aligner_name == "prosplign") {
            tblastn_align(genome_asn, proteins_asn, genome_blastdb, task_params.get('tblastn_align', [:]))
            tblastn_aligns = tblastn_align.out.blast_asn.collect()
            prosplign_prepare(genome_asn, proteins_asn, tblastn_aligns, gencoll_asn, max_intron, task_params.get('prosplign_prepare', [:]))
            prosplign_wnode(genome_asn, proteins_asn, prosplign_prepare.out.compartments_asn, max_intron, task_params.get('prosplign_wnode', [:]))
            alignments_to_use = prosplign_wnode.out.prosplign_asn.collect()
        } else {
            error "Not implemented"
        }

        best_aligned_prot(genome_asn, proteins_asn, alignments_to_use, gencoll_asn, task_params.get('best_aligned_prot', [:]))
        align_filter_sa(genome_asn, proteins_asn, best_aligned_prot.out.asn_file, task_params.get('align_filter_sa', [:]))
        align_sort_sa(genome_asn, proteins_asn,align_filter_sa.out.filtered_file, task_params.get('align_sort_sa', [:]))
    emit:
        protein_alignments = align_sort_sa.out
}
