#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow align_filter_sa {
    take:
        genome_asn_file     //path: genome asn file
        proteins_asn_file   //path: protein asn file
        alignment_asn_file  //path: asn file from best_aligned_prot
        parameters          // Map : extra parameter and parameter update
    main:
        default_params = ""
        effective_params = merge_params(default_params, parameters, 'align_filter')
        run_align_filter_sa(genome_asn_file, proteins_asn_file, alignment_asn_file, effective_params)

    emit:
        filtered_file = run_align_filter_sa.out.filtered_file
        non_match_file = run_align_filter_sa.out.non_match_file
        report_file = run_align_filter_sa.out.report_file
}


process run_align_filter_sa {
    input:
        path genome, stageAs: 'LDS_Index/genome.asnt'
        path proteins,  stageAs: 'LDS_Index/proteins.asnt'
        path asn_file
        val parameters
    output:
        path 'output/align.asn', emit: 'filtered_file'
        path 'output/align-nomatch.asn', emit: 'non_match_file'
        path 'output/report.txt', emit: 'report_file'
        
    script:
    """
    mkdir -p output
    lds2_indexer -source LDS_Index
    align_filter ${parameters}  -lds2 LDS_Index/lds2.db  -nogenbank -input $asn_file -output output/align.asn -non-match-output output/align-nomatch.asn -report-output output/report.txt 
    """
    stub:
    """
    mkdir -p output
    touch output/align.asn
    touch output/align-nomatch.asn
    touch output/report.txt
    """
}
