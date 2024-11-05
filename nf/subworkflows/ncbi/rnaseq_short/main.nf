#!/usr/bin/env nextflow
// rnaseq short  EGAPx execution
// route data to tasks

nextflow.enable.dsl=2

params.import_prefix = "../../../../nf/subworkflows/ncbi/" // redirected during testing

include { sra_query } from "./${params.import_prefix}rnaseq_short/sra_qry/main"
include { fetch_sra_fasta } from "./${params.import_prefix}rnaseq_short/fetch_sra_fasta/main"
include { star_index } from "./${params.import_prefix}rnaseq_short/star_index/main"
include { star_wnode as star } from "./${params.import_prefix}rnaseq_short/star_wnode/main"
include { bam_strandedness } from "./${params.import_prefix}rnaseq_short/bam_strandedness/main"
include { bam_bin_and_sort } from "./${params.import_prefix}rnaseq_short/bam_bin_and_sort/main"
include { bam2asn } from "./${params.import_prefix}rnaseq_short/convert_from_bam/main"
include { rnaseq_collapse } from "./${params.import_prefix}rnaseq_short/rnaseq_collapse/main"


params.intermediate = false


workflow rnaseq_short_plane {
    take:
        genome_asn
        scaffolds
        unpacked_genome_fasta

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
        max_intron      // max intron length
        task_params     // task parameters for every task
    main:
        // Satisfy quirks of Nextflow compiler
        def reads_query1 = reads_query
        def reads_ids1 = reads_ids
        def ch_reads = Channel.fromList(reads)
        // Conditional code on SRA reads source
        if (reads_query || reads_ids || reads) {
            def index = star_index(unpacked_genome_fasta, task_params.get('star_index', [:]))
            def ch_align, ch_align_index, sra_metadata, sra_run_list
            if (reads_query || reads_ids) {
                def query = reads_query1 ? reads_query1 : reads_ids1.join("[Accession] OR ") + "[Accession]"
                (sra_metadata, sra_run_list) = sra_query(query, task_params.get('sra_qry', [:]))
                def reads_fasta_pairs = fetch_sra_fasta(sra_run_list, task_params.get('fetch_sra_fasta', [:]))
                (ch_align, ch_align_index) = star(scaffolds, reads_fasta_pairs, genome_asn, index, max_intron, task_params.get('star_wnode', [:]))
            } else if (ch_reads) {
                sra_metadata = reads_metadata
                (ch_align, ch_align_index) = star(scaffolds, ch_reads, genome_asn, index, max_intron, task_params.get('star_wnode', [:]))
            } else {
                sra_metadata = reads_metadata
                (ch_align, ch_align_index) = star(scaffolds, reads, genome_asn, index, max_intron, task_params.get('star_wnode', [:]))
            }
            //

            bam_strandedness(ch_align.collect(), sra_metadata, task_params.get('bam_strandedness', [:]))
            def strandedness = bam_strandedness.out.strandedness
            
            // Run bam_bin_and_sort
            bam_bin_and_sort(ch_align, ch_align_index, unpacked_genome_fasta, organelles, task_params.get('bam_bin_and_sort', [:]))
            def bam_bins = bam_bin_and_sort.out.sorted

            // Run BAM2ASN
            bam2asn(bam_bins, strandedness, genome_asn, task_params.get('convert_from_bam', [:]))
            def asn_align = bam2asn.out.align.collect()
            def keylist = bam2asn.out.keylist.collect()

            rnaseq_collapse(genome_asn, keylist, asn_align, sra_metadata, 10, task_params.get('rnaseq_collapse', [:]))
        }
        
    emit:
        rnaseq_alignments = rnaseq_collapse.out.alignments 
}
