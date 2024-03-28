#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'

sra_api_url = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo"


workflow sra_query  {
    take:
        // The SRA query should be present either in parameters map in sra_query section or in query
        // The parameter 'query' overrides the one set in 'parameters'
        query          // String : query
        parameters     // Map : extra parameter and parameter update
    main:
        String params = merge_params("", parameters, "sra_query")
        if (query) {
            params = merge_params(params, ["sra_query" : "-query '${query}'"], "sra_query")
        }
        run_sra_query(params)

    emit:
        sra_metadata = run_sra_query.out.sra_metadata
        sra_run_list = run_sra_query.out.sra_run_list
}


process run_sra_query {
    input:
        val parameters
    output:
        path './sra_metadata.dat', emit: 'sra_metadata'
        path './sra_run_accessions.ids', emit: 'sra_run_list'

    script:
    """
    rm -f ./logfile  ./sra_query.ini
    echo -en "[sra]\\napi_url=" > ./sra_query.ini
    echo "$sra_api_url" >> ./sra_query.ini

    sra_query -conffile ./sra_query.ini -logfile ./logfile  -meta ./sra_metadata.dat -output ./sra_run_accessions.ids $parameters 
    rm -f ./logfile  ./sra_query.ini
    """
}