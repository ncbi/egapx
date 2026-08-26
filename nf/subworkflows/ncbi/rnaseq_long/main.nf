#!/usr/bin/env nextflow
// rnaseq long EGAPx execution
// route data to tasks

nextflow.enable.dsl=2

params.import_prefix = "../../../../nf/subworkflows/ncbi/" // redirected during testing

// this process doesnt exist under test
include { rename_fasta_ids } from "./../setup/main"
include { fetch_sra_fasta } from "./${params.import_prefix}rnaseq_short/fetch_sra_fasta/main"
include { minimap2_index } from "./${params.import_prefix}rnaseq_long/minimap2_index/main"
include { minimap2_fasta } from "./${params.import_prefix}rnaseq_long/minimap2_wnode/main"
include { filter_est_align } from "./${params.import_prefix}rnaseq_long/filter_est_align/main"
include { long_read_align_report } from "./${params.import_prefix}rnaseq_long/long_read_align_report/main"
include { gp_register_stats } from '../shared/gp_register_stats/main.nf'

include { checkpoint_save; checkpoint_load_channels; CHECKPOINT_LOAD_JSON} from "./../../../../nf/lib/checkpoint_both"


workflow rnaseq_long_plane_cached {
    take:
        genome_fasta
        gencoll
        // Alternative groups of parameters, one of them should be set
        // reads_ids - list of SRA IDs
        // reads, reads_metadata - path to reads accompanied by metadata
        reads_ids       // list of SRA IDs
        reads           // reads files formatted as fromFilePairs - (list of) tuple(s) [ run_name, [ first_read_file, second_read_file ]]
        max_intron      // max intron length
        long_reads_metadata  // path to reads metadata 13 tab-delimited fields, 1-st - SRA ID, 3-rd paired or unpaired, everything else - not used, but must be present
        task_params     // task parameters for every task
    main:
        def checkpoint_dir = task_params.get('checkpoints_dir', [])
        def checkpoint_name = 'rnaseq_long_plane'
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
            rnaseq_long_plane(genome_fasta, gencoll, reads_ids, reads, max_intron, long_reads_metadata,task_params)
            out_obj = rnaseq_long_plane.out
            checkpoint_save([rnaseq_long_plane.out], checkpoint_args)
        }
    emit:
        alignments = out_obj.alignments
        minimap2_stats = out_obj.minimap2_stats
        align_report = out_obj.align_report
        filter_est_align_stats = out_obj.filter_est_align_stats
}

workflow rnaseq_long_plane {
    take:
        genome_fasta
        gencoll
        // Alternative groups of parameters, one of them should be set
        // reads_ids - list of SRA IDs
        // reads, reads_metadata - path to reads accompanied by metadata
        reads_ids       // list of SRA IDs (ALWAYS a plain list; used as the mode selector)
        reads           // reads as fromFilePairs output: tuple(s) [ run_name, [ read1, read2 ] ]
                        // NOTE: may arrive as a plain List OR as a channel, depending on caller
        max_intron      // max intron length
        long_reads_metadata  // path to reads metadata, 13 tab-delimited fields
        task_params     // task parameters for every task
    main:
        def genome_index = minimap2_index(genome_fasta, task_params.get('minimap2_index', [:]))

        // reads_ids is always a plain list, so this is a safe synchronous emptiness test
        // and serves as the primary mode selector.
        def sra_run_list = reads_ids ? Channel.fromList(reads_ids) : null

        // Tracks whether rename_fasta_ids was actually INVOKED. Gate all .out access on
        // this, never on the reads value (a channel can't be tested synchronously).
        def used_rename = false

        def minimap_wnode_params = task_params.get('minimap2_wnode', [:])
        def fetch_sra_params = task_params.get('fetch_sra_fasta', [:])
        // Pass split from minimap2_wnode to fetch_sra_fasta so splitting happens at download time
        if (minimap_wnode_params.get('split', '')) {
            fetch_sra_params = fetch_sra_params + [split: minimap_wnode_params.get('split')]
        }

        // Does the caller intend the reads path? reads_ids takes precedence.
        // For a List we can test emptiness now; a channel is treated as "present" and
        // will naturally yield nothing downstream if it emits nothing.
        def use_reads = !reads_ids && ( (reads instanceof List) ? !reads.isEmpty() : (reads != null) )

        if (reads_ids) {
            def reads_fasta_pairs = fetch_sra_fasta(sra_run_list, fetch_sra_params)
            // Splitting already happened in fetch_sra_fasta; strip split for the wnode stage.
            def minimap_wnode_params_nosplit = minimap_wnode_params.findAll { k, v -> k != 'split' }
            minimap2_fasta(genome_fasta, genome_index, gencoll, reads_fasta_pairs, max_intron, minimap_wnode_params_nosplit)

        } else if (use_reads) {
            // Normalize reads (List or channel) into a channel of pairs, then assign a
            // REPRODUCIBLE 1-based counter by sorting on the pair name. Sorting makes the
            // counter independent of run-to-run emission order, so a given pair always gets
            // the same srr_id (important for -resume and reproducibility).
            def reads_ch = (reads instanceof List) ? Channel.fromList(reads) : reads

            def indexed = reads_ch
                .toList()                                   // gather all pairs (deferred)
                .flatMap { pairs ->
                    pairs
                        .sort { a, b -> a[0] <=> b[0] }     // deterministic order by run_name
                        .indexed(1)                         // 1-based counter
                        .collect { i, pair -> tuple(pair, i) }
                }

            // Two positional inputs to rename_fasta_ids, derived from one ordered source so
            // they stay aligned emission-for-emission.
            def pairs_ch = indexed.map { pair, i -> pair }   // [ run_name, [files] ]
            def index_ch = indexed.map { pair, i -> i }      // numeric srr_id

            rename_fasta_ids(pairs_ch, index_ch)
            used_rename = true

            minimap2_fasta(genome_fasta, genome_index, gencoll, rename_fasta_ids.out.fasta_pair_list, max_intron, minimap_wnode_params)

        } else {
            minimap2_fasta(genome_fasta, genome_index, gencoll, reads, max_intron, minimap_wnode_params)
        }

        // Merge renamed run metadata with original metadata for long_read_align_report.
        // Gated on used_rename (actual invocation), NOT on the reads value.
        def effective_metadata = long_reads_metadata
        def id_map = Channel.of([])
        if (used_rename) {
            def renamed_metadata = rename_fasta_ids.out.metadata.collectFile(name: 'renamed_runs_metadata.dat', newLine: false)
            effective_metadata = merge_sra_metadata(long_reads_metadata, renamed_metadata)
            id_map = rename_fasta_ids.out.id_map.collectFile(name: 'run_id_map.tsv', newLine: false)
        }

        filter_est_align(minimap2_fasta.out.alignments, gencoll, task_params.get('filter_est_align', [:]))
        long_read_align_report(gencoll, minimap2_fasta.out.minimap2_stats, effective_metadata, filter_est_align.out.stats, id_map, task_params.get('long_read_align_report', [:]))
    emit:
        alignments = filter_est_align.out.alignments
        minimap2_stats = minimap2_fasta.out.minimap2_stats
        align_report = long_read_align_report.out.align_report
        filter_est_align_stats = filter_est_align.out.stats
}


process merge_sra_metadata {
    label 'single_cpu'
    label 'small_mem'
    input:
        path original_metadata
        path renamed_metadata
    output:
        path 'merged_sra_metadata.dat'
    script:
    """
    cat ${original_metadata} > merged_sra_metadata.dat
    cat ${renamed_metadata} >> merged_sra_metadata.dat
    """
    stub:
    """
    cat ${original_metadata} > merged_sra_metadata.dat
    """
}
