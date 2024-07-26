#!/usr/bin/env nextflow
// main nextflow script for EGAPx execution
// route data to subworkflows

nextflow.enable.dsl=2

include { miniprot } from './miniprot/main'
include { align_filter_sa } from './align_filter_sa/main'
include { best_aligned_prot } from './best_aligned_prot/main'
include { paf2asn } from './paf2asn/main'
include { run_align_sort} from '../default/align_sort_sa/main'

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
        run_align_sort(genome_asn, proteins_asn,align_filter_sa.out.filtered_file,
            "-k subject,subject_start,-subject_end,subject_strand,query,query_start,-query_end,query_strand,-num_ident,gap_count" )
    emit:
        protein_alignments = run_align_sort.out
}
