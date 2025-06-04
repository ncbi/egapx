#!/usr/bin/env nextflow
// rnaseq long EGAPx execution
// route data to tasks

nextflow.enable.dsl=2

params.import_prefix = "../../../../nf/subworkflows/ncbi/" // redirected during testing

include { sra_query } from "./${params.import_prefix}rnaseq_short/sra_qry/main"
include { fetch_sra_fasta } from "./${params.import_prefix}rnaseq_short/fetch_sra_fasta/main"
include { minimap2_index } from "./${params.import_prefix}rnaseq_long/minimap2_index/main"
include { minimap2 } from "./${params.import_prefix}rnaseq_long/minimap2_wnode/main"
include { filter_est_align } from "./${params.import_prefix}rnaseq_long/filter_est_align/main"


workflow rnaseq_long_plane {
    take:
        genome_fasta
        gencoll
        // Alternative groups of parameters, one of them should be set
        // reads_query - SRA query in the form accepted by NCBI
        // reads_ids - list of SRA IDs
        // reads, reads_metadata - path to reads accompanied by metadata
        reads_query     // SRA query
        reads_ids       // list of SRA IDs
        reads           // path to reads
        max_intron      // max intron length
        task_params     // task parameters for every task
    main:
        def genome_index = minimap2_index(genome_fasta, task_params.get('minimap2_index', [:]))
        // Satisfy quirks of Nextflow compiler
        def reads_query1 = reads_query
        def reads_ids1 = reads_ids
        //
        def ch_reads = Channel.fromList(reads)
        def sra_metadata, sra_run_list
        if (reads_query || reads_ids) {
            def query = reads_query1 ? reads_query1 : reads_ids1.join("[Accession] OR ") + "[Accession]"
            (sra_metadata, sra_run_list) = sra_query(query, task_params.get('sra_qry', [:]))
            def reads_fasta_pairs = fetch_sra_fasta(sra_run_list, task_params.get('fetch_sra_fasta', [:]))
            minimap2(genome_fasta, genome_index, gencoll, reads_fasta_pairs, max_intron, task_params.get('minimap2_wnode', [:]))
        } else if (ch_reads) {
            minimap2(genome_fasta, genome_index, gencoll, ch_reads, max_intron, task_params.get('minimap2_wnode', [:]))
        } else {
            minimap2(genome_fasta, genome_index, gencoll, reads, max_intron, task_params.get('minimap2_wnode ', [:]))
        }
        filter_est_align(minimap2.out.alignments, task_params.get('filter_est_align', [:]))
    emit:
        alignments = filter_est_align.out.alignments
}
