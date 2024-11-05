#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'

// extract_products -input-manifest models.mft -it -ifmt seq-annot -rna-ids out/rna.ids -prot-ids out/prot.ids

workflow extract_products_from_models {
    take:
        models     //path: models
        parameters          // Map : extra parameter and parameter update
    main:
        default_params = ""
        effective_params = merge_params(default_params, parameters, 'extract_products')
        run_extract_products_from_models(models, effective_params)

    emit:
        rna_ids = run_extract_products_from_models.out.rna_ids
        prot_ids = run_extract_products_from_models.out.prot_ids
        all = run_extract_products_from_models.out.all
}


process run_extract_products_from_models {
    input:
        path models
        val parameters
    output:
        path ('output/rna.ids'), emit: 'rna_ids'
        path ('output/prot.ids'), emit: 'prot_ids'
        path ('output/*'), emit: 'all'
    script:
    """
    mkdir -p output
    echo "${models.join('\n')}"  > models.mft
    extract_products -input-manifest models.mft -it -ifmt seq-annot -rna-ids output/rna.ids -prot-ids output/prot.ids

    """
    stub:
    """
    mkdir -p output
    touch output/rna.ids
    touch output/prot.ids
    """
}
