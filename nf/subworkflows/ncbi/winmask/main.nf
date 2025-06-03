#!/usr/bin/env nextflow


nextflow.enable.dsl=2

params.import_prefix = "../../../../nf/subworkflows/ncbi/" // redirected during testing

include { winmask_wnode } from "./${params.import_prefix}winmask/winmask_wnode/main"
include { make_winmask_stats } from "./${params.import_prefix}winmask/make_winmask_stats/main"
include { prepare_masks } from "./${params.import_prefix}winmask/prepare_masks/main"
include { gc_makeblastdb } from "./${params.import_prefix}winmask/gc_makeblastdb/main"


params.intermediate = false

workflow winmask_plane {
    take:
        genome_asnb
        seqids
        gencoll_asn
        cmsearch_annot
        dustmask_data
        task_params     // task parameters for every task
    main:
        winmask_stats = make_winmask_stats(genome_asnb, seqids, gencoll_asn, task_params.get('make_winmask_stats', [:])).collect()
        mask = winmask_wnode(genome_asnb, seqids, gencoll_asn, winmask_stats.flatten(), task_params.get('winmask_wnode', [:])).collect()
        prepare_masks(cmsearch_annot, dustmask_data, mask, task_params.get('prepare_masks', [:]))
        gc_makeblastdb(genome_asnb, seqids, gencoll_asn, prepare_masks.out.default_softmask_makeblastdb, task_params.get('gc_makeblastdb', [:]))
    emit:
        blastdb = gc_makeblastdb.out.blastdb
        softmask = prepare_masks.out.annot_softmask_gnomon
}
