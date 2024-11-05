#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'

///netmnt/vast01/gpi/regr/GPIPE_REGR1/system/2024-05-16.prod.build26027/bin/compare_annots 
//  -q_scope_type ID 
//  -q_scope_args '' 
//  -q_ids ./gc_get_molecules.8280202/out/no_organelle.gi 
//  -alns ./locus_track_compare.8281472/inp/prev_build_align.mft 
//  -s_scope_type annots-mft-lds 
//  -s_scope_args ./locus_track_compare.8281472/inp/curr_annotation.mft 
//  -s_ids ./passthrough_gcaccess.8279042/out/no_organelle.gi 
//  -o_asn ./locus_track_compare.8281472/out/curr_prev.antcmp.asn 
//  -o_tab ./locus_track_compare.8281472/out/curr_prev.antcmp.tab

///netmnt/vast01/gpi/regr/GPIPE_REGR1/system/2024-05-16.prod.build26027/bin/locus_track 
//   -o_loci      ./locus_track_compare.8281472/out/locus_track.rpt 
//   -o_conflicts ./locus_track_compare.8281472/out/conflicts.rpt 
//   -o_evidence  ./locus_track_compare.8281472/out/evidence.rpt 
//   -annotset (mft of annot_builder ACCEPT files)  curr_annot_set.mft 
//   -gc  (annot assembly, out gc_create) gencoll.asn 
//   -pb_alns (no pb, empty)  prev_build_align.mft 
//   -lxr ./locus_track_compare.8281472/out/lxr_tracking_data.txt 
////     LXR is a File_in to locus_tracki, its from a database call that happens in locus_track's action node, so its saved to the out dir
//   -track_Assemblies var/intra_assm.mft 
//   -track_PrevBuild  (no pb) var/curr_prev.mft

workflow locus_track {
    take:
        annotation 
        gencoll_asn
        lxr_data
        parameters  // Map : extra parameter and parameter update
    main:
        default_params = ""
        effective_params = merge_params(default_params, parameters, 'locus_track')
        run_locus_track(annotation, gencoll_asn, lxr_data, default_params)
    emit:
        track_rpt = run_locus_track.out.track_rpt
}



process run_locus_track {
    input:
        path annotation 
        path gencoll_asn
        path lxr_data //, stageAs: 'lxr_tracking_data.txt'
        val parameters
    output:
        path ('output/locus_track.rpt'), emit: track_rpt
    script:
    """
    mkdir -p output
    ##mkdir -p ./asncache/
    ##prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ....  -oseq-ids /dev/null -split-sequences
    
    echo "${annotation.join('\n')}" > annotation.mft

    locus_track  -nogenbank   -annotset annotation.mft -gc $gencoll_asn -lxr $lxr_data  -o_loci ./output/locus_track.rpt -o_conflicts ./output/conflicts.rpt -o_evidence ./output/evidence.rpt 
    """
    stub:
    """
    mkdir -p output
    touch output/locus_track.rpt
    touch output/lxr_tracking_data.txt
    """
}

