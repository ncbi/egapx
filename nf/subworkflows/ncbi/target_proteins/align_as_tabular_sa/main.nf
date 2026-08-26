#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow align_as_tabular_sa {
    take:
        genome_asnb     //path: genome file
        proteins_asnb   //path: protein file
        asn_file       //path: alignment asn file from paf2asn
        parameters     // Map : extra parameter and parameter update
    main:
        String param = merge_params("", parameters, "align_format")
        run_align_as_tabular_sa(genome_asnb, proteins_asnb, asn_file, param)
    emit:
        align_tab = run_align_as_tabular_sa.out.align_tab
}


process run_align_as_tabular_sa{
    label 'single_cpu'
    label 'small_mem'
    input:
        path genome, stageAs: 'indexed/genome.asnt'
        path proteins,  stageAs: 'indexed/proteins.asnt'
        path asn_file
        val parameters
    output:
        path "output/align.tab", emit: 'align_tab'
    script:
    """
    mkdir -p output
    mkdir -p tmp
    lds2_indexer -source indexed -db tmp/lds_index
    echo "${asn_file.join('\n')}" > align.mft
    align_format -nogenbank -lds2 tmp/lds_index -input-manifest align.mft -o output/align.tab ${parameters}
    rm -rf tmp
    """
    stub:
    """
    mkdir -p output
    touch output/align.tab
    """
}
