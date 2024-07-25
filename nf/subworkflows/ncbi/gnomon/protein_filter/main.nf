#!/usr/bin/env nextflow
nextflow.enable.dsl=2



/*
 *   align_filter -filter 'pct_coverage >= 50' -nogenbank 
 *      | align_sort -ifmt seq-align-set -k query,-bit_score,slen,-align_length -group 1 -top 1 -nogenbank
 */


include { merge_params; to_map; shellSplit } from '../../utilities'


workflow best_protein_hits {
    take:
        gnomon_prot_asn
        swiss_prot_asn
        prot_alignments
        parameters  // Map : extra parameter and parameter update
    main:
        String align_filter_params = merge_params(' -ifmt seq-align-set -filter \'pct_coverage >= 50\' -nogenbank', parameters, 'align_filter')
        String align_sort_params = merge_params(' -ifmt seq-align-set -k query,-bit_score,slen,-align_length -group 1 -top 1 -nogenbank', parameters, 'align_sort')

        run_protein_filter_replacement(gnomon_prot_asn, swiss_prot_asn, prot_alignments, align_filter_params, align_sort_params)

    emit:
        alignments = run_protein_filter_replacement.out
}


process run_protein_filter_replacement {
    input:
        path gnomon_prot_asn, stageAs: 'indexed/*'
        path swiss_prot_asn, stageAs: 'indexed/*'
        path input_prot_alignments, stageAs: "input_alignments.asnb"
        val align_filter_params
        val align_sort_params
    output:
        path "output/*"
    script:
    """
    mkdir -p ./asncache/
    prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ${gnomon_prot_asn} -oseq-ids /dev/null -split-sequences
    prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ${swiss_prot_asn} -oseq-ids /dev/null -split-sequences

    mkdir -p ./output

    align_filter $align_filter_params -asn-cache ./asncache  -i ./input_alignments.asnb -o - | align_sort -i - $align_sort_params -asn-cache ./asncache  -o - | align_pack -ifmt seq-align -i - -ofmt seq-align-set -o ./output/best_protein_hits.asnb 
    ##align_filter $align_filter_params -asn-cache ./asncache  -i ./input_alignments.asnb -o ./t1.asnb 
    ##align_sort -i ./t1.asnb $align_sort_params -asn-cache ./asncache  -o ./t2.asnb     
    ##align_pack -ifmt seq-align -i ./t2.asnb -ofmt seq-align-set -o ./t3.asnb
    ##cp ./t3.asnb ./output/best_protein_hits.asnb

    """

    stub:
    """
    mkdir -p output
    touch ./output/best_protein_hits.asnb
    """
}


