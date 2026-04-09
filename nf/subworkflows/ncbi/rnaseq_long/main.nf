#!/usr/bin/env nextflow
// rnaseq long EGAPx execution
// route data to tasks

nextflow.enable.dsl=2

params.import_prefix = "../../../../nf/subworkflows/ncbi/" // redirected during testing

// this process doesnt exist under test
//include { rename_fasta_ids } from "./${params.import_prefix}/setup/main"
include { rename_fasta_ids } from "./../setup/main"
include { fetch_sra_fasta } from "./${params.import_prefix}rnaseq_short/fetch_sra_fasta/main"
include { minimap2_index } from "./${params.import_prefix}rnaseq_long/minimap2_index/main"
include { minimap2; minimap2_fasta } from "./${params.import_prefix}rnaseq_long/minimap2_wnode/main"
include { filter_est_align } from "./${params.import_prefix}rnaseq_long/filter_est_align/main"


workflow rnaseq_long_plane {
    take:
        genome_fasta
        gencoll
        // Alternative groups of parameters, one of them should be set
        // reads_ids - list of SRA IDs
        // reads, reads_metadata - path to reads accompanied by metadata
        reads_ids       // list of SRA IDs
        reads           // reads files formatted as fromFilePairs - (list of) tuple(s) [ run_name, [ first_read_file, second_read_file ]]
        max_intron      // max intron length
        task_params     // task parameters for every task
    main:
        def genome_index = minimap2_index(genome_fasta, task_params.get('minimap2_index', [:]))
        // Satisfy quirks of Nextflow compiler
        def reads_ids1 = reads_ids
        def sra_run_list = reads_ids ? Channel.fromList(reads_ids1) : null
        def ch_reads = null
        def reads1 = reads
        def reads_indices = null
        try {
            ch_reads = Channel.fromList(reads)
            reads_indices = Channel.from(1..reads.size())
        } catch( Exception e ) {
            ch_reads = null
        }
        def minimap_wnode_params = task_params.get('minimap2_wnode', [:])
        def fetch_sra_params = task_params.get('fetch_sra_fasta', [:])
        // Pass split parameter from minimap2_wnode to fetch_sra_fasta so splitting happens at download time
        if (minimap_wnode_params.get('split', '')) {
            fetch_sra_params = fetch_sra_params + [split: minimap_wnode_params.get('split')]
        }

        if (reads_ids) {
            def reads_fasta_pairs = fetch_sra_fasta(sra_run_list, fetch_sra_params)
            // Remove split from minimap_wnode_params since splitting already happened in fetch_sra_fasta
            def minimap_wnode_params_nosplit = minimap_wnode_params.findAll { k, v -> k != 'split' }
            alignments = minimap2_fasta(genome_fasta, genome_index, gencoll, reads_fasta_pairs, max_intron, minimap_wnode_params_nosplit)
        } else if (ch_reads) {
            def renamed_reads = rename_fasta_ids(ch_reads, reads_indices)
            // Use minimap2_fasta for explicit reads too (supports optional splitting)
            alignments = minimap2_fasta(genome_fasta, genome_index, gencoll, renamed_reads, max_intron, minimap_wnode_params)
        } else {
            alignments = minimap2_fasta(genome_fasta, genome_index, gencoll, reads, max_intron, minimap_wnode_params)
        }
        filter_est_align(alignments, task_params.get('filter_est_align', [:]))
    emit:
        alignments = filter_est_align.out.alignments
}
