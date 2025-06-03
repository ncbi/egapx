#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow minimap2_index {
    take:
        genome
        parameters  // Map : extra parameter and parameter update
    main:
        String index_params =  merge_params("-k 15 -w 5", parameters, 'minimap2_index')
        run_minimap2_index(genome, index_params)
    emit:
        index = run_minimap2_index.out
}


process run_minimap2_index {
    input:
        path genome
        val index_params
    output:
        path 'index.mmi'
    script:
    """
    minimap2 -d index.mmi $index_params $genome
    """
    stub:
    """
    touch index.mmi
    """
}