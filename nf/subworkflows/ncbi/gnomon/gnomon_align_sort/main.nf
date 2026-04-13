#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { merge_params } from '../../utilities'
include { run_align_sort }  from '../../default/align_sort_sa/main.nf'

workflow gnomon_align_sort {
    take:
        alignments
        parameters  // Map : extra parameter and parameter updates
    main:
        String input_sorting = parameters.get('input_aligns_sort', '')
        def sort_aligns = alignments
        if (!input_sorting.contains("presorted")) {
            String lcl_params
            if (input_sorting.contains("merge_only")) {
                lcl_params = "-merge"
            }
            lcl_params += " -nogenbank -ifmt seq-align -compression none -k 'subject,subject_start,-subject_end' -rnaseq-uniq -strip-alignment "
            String align_sort_params = merge_params(lcl_params, parameters, 'align_sort')
            sort_aligns = run_align_sort([], [], alignments, align_sort_params).collect()
        }
    emit:
        sorted_alignments = sort_aligns
}
