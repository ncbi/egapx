#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow prot_gnomon_prepare {
    take:
        models
        parameters  // Map : extra parameter and parameter update
    main:
        String prot_gnomon_prepare_params =  merge_params('', parameters, 'prot_gnomon_prepare')

        prot_gnomon_prepare_p(models, prot_gnomon_prepare_params)
    emit:
        outputs = prot_gnomon_prepare_p.out.outputs
        LDS2 = prot_gnomon_prepare_p.out.LDS2
        prot_ids = prot_gnomon_prepare_p.out.prot_ids
        nuc_ids = prot_gnomon_prepare_p.out.nuc_ids
}


process prot_gnomon_prepare_p {
    input:
        path models
        val params
    output:
        path "*", emit: 'outputs'
        path "LDS2", emit: 'LDS2'
        path "prot_ids.seq_id", emit: 'prot_ids'
        path "nuc_ids.seq_id", emit: 'nuc_ids'
    script:
    """
    echo ${models} |tr ' ' '\\n' > models.mft
    prot_gnomon_prepare ${params} -input-manifest models.mft -olds2 LDS2 -oprot-ids prot_ids.seq_id -onuc-ids nuc_ids.seq_id
    """
    stub:
    """
    touch LDS2
    touch prot_ids.seq_id
    touch nuc_ids.seq_id
    """
}

