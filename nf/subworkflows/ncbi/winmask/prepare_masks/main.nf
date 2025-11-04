#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params; convert_mask as convert_rfam_mask; 
        convert_mask as convert_softmask; combine_blast_db as combine_blast_db_softmask;
        combine_blast_db as combine_blast_db_alt_softmask } from '../../utilities'


workflow prepare_masks {
    take:
        cmsearch_annot
        //contam_hardmask
        //contam_softmask
        //rmask_data
        dustmask_data
        winmask_data
        //rrna_mask_data
        parameters
    main:
        rfam_rrna =  gp_annot_format(cmsearch_annot)
        String convert_mask_params = merge_params("", parameters, "convert_mask")
        String combine_blast_db_params = merge_params("", parameters, "combine_blast_db")
        rfam_rrna_masks = convert_rfam_mask(rfam_rrna,  "rfam_rrna_masks", convert_mask_params)
        //String combine_blast_db_hardmask_params = merge_params("", parameters, "combine_blast_db_hardmask")
        String combine_blast_db_softmask_params = merge_params("", parameters, "combine_blast_db_softmask")
        //default_hardmask = combine_blast_db(contam_hardmask, [], [], [], [],[], combine_blast_db_hardmask_params)
        default_softmask = combine_blast_db_softmask(dustmask_data, winmask_data,/* rrna_mask_data, */ "default_softmask", combine_blast_db_softmask_params)
        alternate_softmask = combine_blast_db_alt_softmask(winmask_data, rfam_rrna_masks, "alternate_softmask", combine_blast_db_params)
        //String convert_mask_hardmask_params = merge_params("", parameters, "convert_mask_hardmask")
        String convert_mask_softmask_params = merge_params("", parameters, "convert_mask_softmask")
        //default_hardmask_asnb = convert_mask(default_hardmask, convert_mask_hardmask_params)
        default_softmask_asnb = convert_softmask(default_softmask, "default_softmask_asnb", convert_mask_softmask_params)
    emit:
        //annot_hardmask = default_hardmask_asnb
        annot_softmask_gnomon = default_softmask_asnb
        //default_hardmask = default_hardmask
        default_softmask_makeblastdb = default_softmask
        alternate_softmask = alternate_softmask
}


process gp_annot_format{
    input:
        path cmsearch_annot
    output:
        path "output/rfam_rrna.tab", emit: 'rfam_rrna'
    script:
    """
    mkdir -p output
    echo "${cmsearch_annot}" > rfam.mft
    gp_annot_format -input-manifest rfam.mft -ifmt seq-entry -ofmt tabular -o output/rfam_rrna.tab -nogenbank
    """
    stub:
    """
    mkdir -p output
    touch output/rfam_rrna.tab
    """
}