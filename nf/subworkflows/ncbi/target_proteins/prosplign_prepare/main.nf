#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow prosplign_prepare {
    take:
        genome_asnb     //path: genome file
        proteins_asnb   //path: protein file
        seed_hits       //path: alignment asn file from paf2asn
        gencoll_asn     //path:
        max_intron      //val
        parameters      // Map : extra parameter and parameter update
    main:
        String prosplign_compart_params = merge_params("", parameters, "prosplign_compart")
        String prosplign_compart_filter_params = merge_params("", parameters, "prosplign_compart_filter")
        run_prosplign_prepare(genome_asnb, proteins_asnb, gencoll_asn, seed_hits, max_intron, prosplign_compart_params, prosplign_compart_filter_params)
    emit:
        compartments_asn = run_prosplign_prepare.out.compartments_asn
}


process run_prosplign_prepare{
    label 'single_cpu'
    label 'med_mem'
    input:
        path genome_asnb
        path proteins_asn
        path gencoll_asn
        path  seed_hits
        val max_intron
        val prosplign_compart_params
        val prosplign_compart_filter_params
    output:
        path "output/compartments.asn", emit: 'compartments_asn'
    script:
    """
    mkdir -p output
    mkdir -p tmp/asncache
    echo "${seed_hits.join('\n')}" > seed_hits.mft
    auto_prime_cache.py -cache tmp/asncache -i ${genome_asnb} -oseq-ids spids ##-split-sequences
    auto_prime_cache.py -cache tmp/asncache -i ${proteins_asn} -oseq-ids spids2 -split-sequences
    prosplign_compart -asn-cache tmp/asncache -nogenbank -max_intron $max_intron -input-manifest seed_hits.mft $prosplign_compart_params -o unfiltered_compartments.asn -ifmt seq-align-set
    prosplign_compart_filter -asn-cache tmp/asncache -nogenbank -i unfiltered_compartments.asn -gc-assembly $gencoll_asn $prosplign_compart_filter_params -o output/compartments.asn
    rm -rf tmp
    """
    stub:
    """
    mkdir -p output
    touch output/compartments.asn
    """
}
