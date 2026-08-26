#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow align_filter_sa {
    take:
        genome_asn_file     //path: genome asn file
        proteins_asn_file   //path: protein asn file
        alignment_asn_file  //path: asn file from best_aligned_prot
        id_files            // path: id files
        parameters          // Map : extra parameter and parameter update
    main:
        default_params = ""
        effective_params = merge_params(default_params, parameters, 'align_filter')
        run_align_filter_sa(genome_asn_file, proteins_asn_file, alignment_asn_file, id_files, effective_params)

    emit:
        filtered_file = run_align_filter_sa.out.filtered_file
        non_match_file = run_align_filter_sa.out.non_match_file
        report_file = run_align_filter_sa.out.report_file
}


process run_align_filter_sa {
    label 'multi_cpu'
    label 'small_mem'
    input:
        path genome, stageAs: 'indexed/genome.asnt'
        path proteins,  stageAs: 'indexed/proteins.asnt'
        path asn_file
        path id_files
        val parameters
    output:
        path 'output/align.asn', emit: 'filtered_file'
        path 'output/align-nomatch.asn', emit: 'non_match_file'
        path 'output/report.txt', emit: 'report_file'
        
    script:
    """
    mkdir -p output
    mkdir -p tmp
    lds2_indexer -source indexed -db tmp/lds_index
    echo "${id_files.join('\n')}" > include_list.mft
    #for align_filter_sa_remove_low_quality, so we need to remove the high_ident_xps.ids from the include list.
    sed -i '/high_ident_xps.ids/d' include_list.mft
    align_filter ${parameters} -query-allowlist include_list.mft -lds2 tmp/lds_index  -nogenbank -input $asn_file -output output/align.asn -non-match-output output/align-nomatch.asn -report-output output/report.txt 
    rm -rf tmp
    """
    stub:
    """
    mkdir -p output
    touch output/align.asn
    touch output/align-nomatch.asn
    touch output/report.txt
    """
}
