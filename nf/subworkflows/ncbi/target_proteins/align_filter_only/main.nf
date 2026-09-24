#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { merge_params } from '../../utilities'
include { run_align_filter_sa as run_align_filter_only }  from '../../default/align_filter_sa/main.nf'

workflow align_filter_only {
    take:
        genome_asnb     // path: genome file
        proteins_asnb   // path: protein file
        alignments        // path: alignment asn file from paf2asn
        id_files        //
        parameters      // Map : extra parameter and parameter update
    main:
        String align_filter_only_params = merge_params("", parameters, 'align_filter_only')
        (filtered, non_match, report) = run_align_filter_only(genome_asnb, proteins_asnb, alignments, id_files, align_filter_only_params)
    emit:
        filtered_file = filtered
        non_match_file = non_match
        report_file = report
}
