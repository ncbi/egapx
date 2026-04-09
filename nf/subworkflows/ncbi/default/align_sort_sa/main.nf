#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { merge_params} from '../../utilities'


workflow align_sort_sa {
    take:
        genome_asn
        proteins_asn 
        alignments  //path: alignment files
        parameters      // Map : extra parameter and parameter update
    main:
        default_params = " -nogenbank "
        effective_params = merge_params(default_params, parameters, 'align_sort')
        run_align_sort(genome_asn, proteins_asn, alignments, effective_params)

    emit:
        sorted_asn_file = run_align_sort.out.sorted_asn_file
}


process run_align_sort {
    label 'multi_cpu'
    label 'large_mem'
    input:
        path genome, stageAs: 'indexed/genome.asnt'
        path proteins, stageAs: 'indexed/proteins.asnt'
        path alignments
        val parameters
    output:
        path ('output/*.asn')  , emit: 'sorted_asn_file'
    script:
    """
    mkdir -p output
    mkdir -p tmp
    lds2_indexer -source indexed -db tmp/lds_index
    echo "${alignments.join('\n')}" > alignments.mft
    align_sort  $parameters -tmp tmp -input-manifest alignments.mft  -o output/sorted_aligns.asn  -lds2 tmp/lds_index
    rm -rf tmp
    """
    stub:
    """
    mkdir -p output
    name=`echo "${alignments.join('\n')}" | md5sum | cut -c-16`
    touch output/\${name}.sorted_aligns.asn
    """
}
