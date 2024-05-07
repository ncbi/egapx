#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'

workflow gnomon_training {
    take:
        genome_asn
        models_file
        parameters  // Map : extra parameter and parameter update
    main:
        default_params = "-b -asn -maxintron 1200000"
        effective_params = merge_params(default_params, parameters, 'gnomon_training')
        run_gnomon_training(genome_asn, models_file, default_params)
    emit:
        hmm_params_file = run_gnomon_training.out.hmm_params_file
}



process run_gnomon_training {
    input:
        path genome_asn, stageAs: 'indexed/*'
        path models_file
        val parameters
    output:
        path ('output/hmm_params.asn'), emit: 'hmm_params_file'

    script:
    """
    mkdir -p output
    lds2_indexer -source indexed -db ./indexed_lds
    gnomon_training ${parameters} -nogenbank -lds2 ./indexed_lds  -input ${models_file} -out output/hmm_params.asn 
    """
    stub:
    """
    mkdir -p output
    touch  output/hmm_params.asn
    """
}
