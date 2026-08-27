#!/usr/bin/env nextflow
// rnaseq short  EGAPx execution
// route data to tasks

nextflow.enable.dsl=2

params.import_prefix = "../../../../nf/subworkflows/ncbi/" // redirected during testing

include { fetch_sra_fasta } from "./${params.import_prefix}rnaseq_short/fetch_sra_fasta/main"
include { star_index } from "./${params.import_prefix}rnaseq_short/star_index/main"
include { star_wnode as star } from "./${params.import_prefix}rnaseq_short/star_wnode/main"
include { bam_strandedness } from "./${params.import_prefix}rnaseq_short/bam_strandedness/main"
include { bam_bin_and_sort } from "./${params.import_prefix}rnaseq_short/bam_bin_and_sort/main"
include { bam2asn } from "./${params.import_prefix}rnaseq_short/convert_from_bam/main"
include { rnaseq_collapse } from "./${params.import_prefix}rnaseq_short/rnaseq_collapse/main"
include { rnaseq_register_stats } from "./${params.import_prefix}rnaseq_short/rnaseq_register_stats/main"
include { rnaseq_align_report } from "./${params.import_prefix}rnaseq_short/rnaseq_align_report/main"

include { checkpoint_save; checkpoint_load_channels; CHECKPOINT_LOAD_JSON} from "./../../../../nf/lib/checkpoint_both"


params.intermediate = false


workflow rnaseq_short_plane_cached {
    take:
        genome_asn
        gencoll_asn
        scaffolds
        unpacked_genome_fasta
        reads_ids       // list of SRA IDs
        reads           // path to reads
        reads_metadata  // path to reads metadata 13 t
        organelles      // path to organelle list
        tax_id          // NCBI tax id of the closest taxon to the genome
        max_intron      // max intron length
        task_params     // task parameters for every task
    main:
        def checkpoint_dir = task_params.get('checkpoints_dir', [])
        def checkpoint_name = 'rnaseq_short_plane'
        def checkpoint_enabled = task_params.get('checkpoints_save', false)
        def checkpoint_args = [name: checkpoint_name, dir: checkpoint_dir, enabled: checkpoint_enabled]
        def ck_file = file("${checkpoint_dir}/${checkpoint_name}.json")
        def out_obj=[:]
        if (ck_file.exists()) {
            CHECKPOINT_LOAD_JSON(checkpoint_name, checkpoint_dir)
            //def loaded = checkpoint_load_channels(CHECKPOINT_LOAD_JSON.out.json_file)
            def loaded = checkpoint_load_channels(checkpoint_args)
            out_obj = loaded
        } else {
            rnaseq_short_plane(genome_asn, gencoll_asn, scaffolds, 
                            unpacked_genome_fasta, reads_ids, reads, reads_metadata,
                            organelles, tax_id, max_intron, task_params)
            out_obj = rnaseq_short_plane.out
            checkpoint_save([rnaseq_short_plane.out], checkpoint_args)
        }
    emit:
        rnaseq_alignments = out_obj.rnaseq_alignments 
        sra_exons = out_obj.sra_exons
        sra_exons_slices = out_obj.sra_exons_slices
        star_bam = out_obj.star_bam
        run_stats = out_obj.run_stats
        align_report = out_obj.align_report
        run_reports = out_obj.run_reports
        star_logs = out_obj.star_logs
}

workflow rnaseq_short_plane {
    take:
        genome_asn
        gencoll_asn
        scaffolds
        unpacked_genome_fasta

        // Alternative groups of parameters, one of them should be set
        // reads_ids - list of SRA IDs
        // reads, reads_metadata - path to reads accompanied by metadata
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
        def reads_ids1 = reads_ids
        def sra_run_list = reads_ids ? Channel.fromList(reads_ids1) : null
        def ch_reads = []
        try {
            ch_reads = Channel.fromList(reads)
        }
        catch( Exception e ) {
           ch_reads = reads
        }

        star_bam_out = []
        star_logs_out = []
        // Conditional code on SRA reads source
        if (reads_ids || reads) {
            def index = star_index(unpacked_genome_fasta, task_params.get('star_index', [:]))
            def ch_align, ch_align_index, ch_star_logs
            if (reads_ids) {
                def reads_fasta_pairs = fetch_sra_fasta(sra_run_list, task_params.get('fetch_sra_fasta', [:]))
                (ch_align, ch_align_index, ch_star_logs) = star(scaffolds, reads_fasta_pairs, genome_asn, index, max_intron, task_params.get('star_wnode', [:]))
            } else if (ch_reads) {
                (ch_align, ch_align_index, ch_star_logs) = star(scaffolds, ch_reads, genome_asn, index, max_intron, task_params.get('star_wnode', [:]))
            } else {
                (ch_align, ch_align_index, ch_star_logs) = star(scaffolds, reads, genome_asn, index, max_intron, task_params.get('star_wnode', [:]))
            }
            //

            bam_strandedness(ch_align.collect(), reads_metadata, task_params.get('bam_strandedness', [:]))
            def strandedness = bam_strandedness.out.strandedness

            // Run bam_bin_and_sort
            bam_bin_and_sort(ch_align, ch_align_index, unpacked_genome_fasta, organelles, task_params.get('bam_bin_and_sort', [:]))
            def bam_bins = bam_bin_and_sort.out.sorted

            // Run BAM2ASN
            bam2asn(bam_bins, strandedness, genome_asn, reads_metadata, task_params.get('convert_from_bam', [:]))
            def asn_align = bam2asn.out.align.collect()
            def keylist = bam2asn.out.keylist.collect()

            rnaseq_collapse(genome_asn, keylist, asn_align, reads_metadata, 10, task_params.get('rnaseq_collapse', [:]))

            // Run rnaseq_register_stats
            def collapsed_aligns = rnaseq_collapse.out.alignments
            def per_run_counts = keylist
            def stranded_runs = bam_strandedness.out.stranded_runs.ifEmpty([])
            def unstranded_runs = bam_strandedness.out.unstranded_runs.ifEmpty([])
            rnaseq_register_stats(gencoll_asn, reads_metadata, collapsed_aligns, per_run_counts, stranded_runs, unstranded_runs, task_params.get('rnaseq_register_stats', [:]))

            // Run rnaseq_align_report
            def run_stats = rnaseq_register_stats.out.outputs
            def run_list = strandedness
                .map { strandedness_file ->
                    strandedness_file
                        .readLines()
                        .findAll { line -> line && !line.startsWith('#') }
                        .collect { line -> line.tokenize('\t')[0] }
                }
                .flatten()
                .unique()
                .collect()

            rnaseq_align_report(gencoll_asn, reads_metadata, per_run_counts, run_stats, run_list, task_params.get('rnaseq_align_report', [:]))
            star_bam_out = ch_align
            star_logs_out = ch_star_logs
        }
        
    emit:
        rnaseq_alignments = rnaseq_collapse.out.alignments 
        sra_exons = rnaseq_collapse.out.exons
        sra_exons_slices = rnaseq_collapse.out.exons_slices
        star_bam = star_bam_out
        run_stats = rnaseq_register_stats.out.outputs
        align_report = rnaseq_align_report.out.align_report
        run_reports = rnaseq_align_report.out.run_reports
        star_logs = star_logs_out
}
