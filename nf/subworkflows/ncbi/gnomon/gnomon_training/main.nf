#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'

workflow gnomon_training {
    take:
        models_file
        parameters  // Map : extra parameter and parameter update
    main:
        default_params = "-b -asn -maxintron 1200000"
        effective_params = merge_params(default_params, parameters, 'gnomon_training')
        run_gnomon_training(models_file, default_params)
    emit:
        hmm_params_file = run_gnomon_training.out.hmm_params_file
}



process run_gnomon_training {
    input:
        path models_file
        val parameters
    output:
        path ('output/hmm_params.asn'), emit: 'hmm_params_file'

    script:
    """
    mkdir -p output
    gnomon_training ${parameters} -input ${models_file} -out output/hmm_params.asn 
    """
}
