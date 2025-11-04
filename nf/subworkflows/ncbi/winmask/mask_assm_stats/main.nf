#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params;} from '../../utilities'


workflow mask_assm_stats {
    take:
        gencoll_asn
        genome_asnb
        winmask
        mask_data
        parameters
    main:
        mask_stats = run_mask_assm_stats(gencoll_asn,genome_asnb, winmask, mask_data)
    emit:
        mask_stats = mask_stats
}


process run_mask_assm_stats{
    input:
        path gencoll_asn
        path genome_asnb
        path winmask
        path mask_data
    output:
        path "output/mask_assm_stats.xml", emit: 'mask_stats'
    script:
    """
    mkdir -p output
    echo $winmask > winmask_data.mft
    echo $mask_data > rmask_data.mft
        mkdir -p ./asncache/
    prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ${genome_asnb} -oseq-ids spids -split-sequences
    mask_assm_stats -assembly $gencoll_asn -nogenbank -rmask-output rmask_data.mft -winmask-output winmask_data.mft -output output/mask_assm_stats.xml -asn-cache ./asncache/
    """
    stub:
    """
    mkdir -p output
    touch output/mask_assm_stats.xml
    """
}