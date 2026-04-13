#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { merge_params } from '../../utilities'

workflow filter_est_align {
    take:
        align_asn
        parameters  // Map : extra parameter and parameter update
    main:
        String align_filter_params = merge_params('-ifmt seq-align -nogenbank', parameters, 'align_filter')
        String align_sort_params = merge_params('-ifmt seq-align -nogenbank', parameters, 'align_sort')
        align_filter(align_asn, align_filter_params)
        align_sort(align_filter.out, align_sort_params)
    emit:
        alignments = align_sort.out
}


process align_filter {
    label 'single_cpu'
    label 'small_mem'
    input:
        path align_asn, stageAs: 'indexed/*'
        val align_filter_params
    output:
        path "output/*"
    script:
    name = align_asn.baseName
    """
    mkdir -p tmp
    mkdir -p output
    lds2_indexer -source ./indexed -db tmp/LDS2_Index
    align_filter $align_filter_params -lds2 tmp/LDS2_Index -i ${align_asn} -o output/${name}.filtered_aln.asnb
    """
    stub:
    name = align_asn.baseName
    """
    mkdir -p output
    touch output/${name}.filtered_aln.asnb
    """
}


process align_sort {
    label 'single_cpu'
    label 'large_mem'
    input:
        path align_asn, stageAs: 'indexed/*'
        val align_sort_params
    output:
        path "output/*"
    script:
    name = align_asn.baseName
    """
    mkdir -p output
    mkdir -p tmp
    lds2_indexer -source ./indexed -db tmp/LDS2_Index
    align_sort -tmp tmp  $align_sort_params -lds2 tmp/LDS2_Index -i ${align_asn} -o output/${name}.sorted_aln.asnb
    rm -rf tmp
    """
    stub:
    name = align_asn.baseName
    println("align_sort parameters: $align_sort_params")
    """
    mkdir -p output
    touch output/${name}.sorted_aln-long-reads.asnb
    """
}
