#!/usr/bin/env nextflow
// main nextflow script for EGAPx execution
// route data to subworkflows

nextflow.enable.dsl=2

params.import_prefix = "../../../../nf/subworkflows/ncbi/" // redirected during testing

include { miniprot } from "./${params.import_prefix}target_proteins/miniprot/main"
include { align_filter_sa } from "./${params.import_prefix}target_proteins/align_filter_sa/main"
include { best_aligned_prot } from "./${params.import_prefix}target_proteins/best_aligned_prot/main"
include { paf2asn } from "./${params.import_prefix}target_proteins/paf2asn/main"
include { align_sort_sa} from "./${params.import_prefix}target_proteins/../default/align_sort_sa/main"

params.intermediate = false


workflow target_proteins_plane {
    take:
        unpacked_genome_fasta
        genome_asn
        gencoll_asn
        unpacked_proteins_fasta
        proteins_asn
        max_intron
        task_params     // task parameters for every task
    main:
        // Protein alignments
        miniprot(unpacked_genome_fasta, unpacked_proteins_fasta, max_intron, task_params.get('miniprot', [:]))
        def miniprot_file = miniprot.out.miniprot_file
        paf2asn(genome_asn, proteins_asn, miniprot_file, task_params.get('paf2asn', [:]))
        def converted_asn = paf2asn.out.asn_file
        best_aligned_prot(genome_asn, proteins_asn, converted_asn.collect(), gencoll_asn, task_params.get('best_aligned_prot', [:]))
        align_filter_sa(genome_asn, proteins_asn, best_aligned_prot.out.asn_file, task_params.get('align_filter_sa', [:]))
        align_sort_sa(genome_asn, proteins_asn,align_filter_sa.out.filtered_file, task_params.get('align_sort_sa', [:]))
    emit:
        protein_alignments = align_sort_sa.out
}
