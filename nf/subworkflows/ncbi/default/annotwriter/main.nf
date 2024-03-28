#!/usr/bin/env nextflow
nextflow.enable.dsl=2


workflow annotwriter {
    take:
        accept_asn_file  // Channel: accept.asn file
        parameters     // Map : extra parameter and parameter update
    main:
        c =  run_annotwriter(accept_asn_file)
    emit:
        annoted_file = c
}


process run_annotwriter {
    input:
        path accept_asn_file
    output:
        path ('output/accept.gff')  , emit: 'annoted_file'

    script:
    """
    mkdir -p output
    annotwriter -i ${accept_asn_file} -nogenbank -format gff3 -o output/accept.gff
    """
}
