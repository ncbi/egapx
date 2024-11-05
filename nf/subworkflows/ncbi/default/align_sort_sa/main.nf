#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { merge_params} from '../../utilities'


workflow align_sort_sa {
    take:
        genome_asn
        proteins_asn 
        alignments  //path: alignment files
        parameters      // Map : extra parameter and parameter update
    main:
        default_params = ""
        effective_params = merge_params(default_params, parameters, 'align_sort')
        run_align_sort( genome_asn, proteins_asn, alignments, effective_params)

    emit:
        sorted_asn_file = run_align_sort.out.sorted_asn_file
}


process run_align_sort {
    input:
        path genome, stageAs: 'LDS_Index/genome.asnt'
        path proteins,  stageAs: 'LDS_Index/proteins.asnt'
        path alignments
        val parameters
    output:
        path ('output/sorted_aligns.asn')  , emit: 'sorted_asn_file'
    script:
    """
    mkdir -p output
    mkdir -p LDS_Index
    lds2_indexer -source LDS_Index
    echo "${alignments.join('\n')}" > alignments.mft
    mkdir -p tmp
    align_sort  $parameters -tmp tmp -input-manifest alignments.mft  -o output/sorted_aligns.asn  -lds2 LDS_Index/lds2.db
    """
    stub:
    """
    mkdir -p output
    touch output/sorted_aligns.asn
    """
}
