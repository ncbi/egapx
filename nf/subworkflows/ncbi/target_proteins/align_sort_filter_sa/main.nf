#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { merge_params } from '../../utilities'
include { run_align_sort }  from '../../default/align_sort_sa/main.nf'

workflow align_sort_filter_sa {
    take:
        genome_asnb     // path: genome file
        proteins_asnb   // path: protein file
        asn_file        // path: alignment asn file from paf2asn
        parameters      // Map : extra parameter and parameter update
    main:
        String align_sort_params = parameters["align_sort"]
        def sort_aligns = run_align_sort(genome_asnb, proteins_asnb, asn_file, align_sort_params).collect()
    emit:
        sorted_asn_file = sort_aligns
}