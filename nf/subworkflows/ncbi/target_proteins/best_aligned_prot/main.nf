#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow best_aligned_prot {
    take:
        genome_asn_file     //path: genome file
        proteins_asn_file   //path: protein file
        alignment_asn_file  //path: alignment asn file from paf2asn
        gencoll_file        //path:
        parameters          // Map : extra parameter and parameter update
    main:
        default_params = ""
        effective_params = merge_params(default_params, parameters, 'best_placement')
        run_best_aligned_prot(genome_asn_file, proteins_asn_file, alignment_asn_file, gencoll_file, effective_params)

    emit:
        asn_file = run_best_aligned_prot.out.asn_file
        report_file = run_best_aligned_prot.out.report_file
        
}


process run_best_aligned_prot {
    input:
        path genome, stageAs: 'LDS_Index/genome.asnt'
        path proteins,  stageAs: 'LDS_Index/proteins.asnt'
        path alignment_asn_file // list of alignment files
        path gencoll_file
        val parameters
    output:
        path 'output/align.asn', emit: 'asn_file'
        path 'output/report.txt', emit: 'report_file'
        
    script:
    """
    mkdir -p output
    lds2_indexer -source LDS_Index
    echo "${alignment_asn_file.join('\n')}" > align.mft
    best_placement ${parameters}  -lds2 LDS_Index/lds2.db  -nogenbank  -gc_path $gencoll_file -in_alns align.mft -out_alns output/align.asn -out_rpt  output/report.txt
    """
    stub:
        print("Best aligned prot input ${alignment_asn_file}")
    """
    mkdir -p output
    touch output/align.asn
    touch output/report.txt
    """
}
