#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow paf2asn {
    take:
        genome_asn_file   //path: genome asn file
        proteins_asn_file //path: protein asn file
        paf_file          //path: paf alignment file from miniprot
        parameters        //Map : extra parameter and parameter update
    main:
        default_params = ""
        effective_params = merge_params(default_params, parameters, 'paf2asn')
        run_paf2asn(genome_asn_file, proteins_asn_file, paf_file, effective_params)

    emit:
        asn_file = run_paf2asn.out.asn_file
}


process run_paf2asn {
    label 'long_job'
    label 'multi_node'
    label 'small_mem'
    input:
        path genome,  stageAs: 'indexed/genome.asnt'
        path proteins,  stageAs: 'indexed/proteins.asnt'
        path paf_file  // list of PAF files to convert
        val parameters
    output:
        path 'output/*.asn', emit: 'asn_file'
    script:
        def asn_name = paf_file.baseName.toString() + ".asn"    // fails when it it multiple files (paf_file)
        //def asn_name = "paf_file.align.asn"    // this could be multiple files (paf_file)
    """
    mkdir -p output
    mkdir -p tmp
    lds2_indexer -source indexed -db tmp/lds2_index
    echo "${paf_file.join('\n')}" > input.mft
    paf2asn ${parameters} -lds2 tmp/lds2_index -nogenbank -input-manifest input.mft |
        align_filter -u -nogenbank | # deduplicate - GP-40550
        cat > output/${asn_name}
    rm -rf tmp
    """
    stub:
        def asn_name = paf_file.baseName.toString() + ".asn"    
    """
    mkdir -p output
    touch output/${asn_name}
    """
}
