#!/usr/bin/env nextflow
nextflow.enable.dsl=2

 
include { merge_params } from '../../utilities'


workflow gnomon_add_filtering_scores {
    take:
        naming_db_secondary_support
        search_set_secondary_support
        gnomon_filtering_scores_file
        gnomon_quality_report
        best_naming_hits
        models
        swiss_prot_asn
        parameters  // Map : extra parameter and parameter update
    main:
        default_params = " -naming-db-support-filter 'pct_coverage >= 70 AND align_length >= 100 AND symmetric_overlap >= 75' -search-set-support-filter '0=1' -cap-and-polya-markup"
        effective_params = merge_params(default_params, parameters, 'gnomon_filter_models')
        filter_models = gnomon_filter_models(naming_db_secondary_support, search_set_secondary_support, gnomon_filtering_scores_file, gnomon_quality_report , best_naming_hits, models, swiss_prot_asn, effective_params)
    emit:
        gnomon_filter_models = filter_models
}

process gnomon_filter_models {
    input:
        path naming_db_secondary_support
        path search_set_secondary_support
        path gnomon_filtering_scores_file 
        path gnomon_quality_report
        path best_naming_hits
        path models
        path swiss_prot_asn
        val params
    output:
        path "scored.models.asn", emit: "gnomon_filter_models"
    script:
    """

    mkdir -p ./asncache/
    prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ${models} -oseq-ids spids -split-sequences
    prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ${swiss_prot_asn} -oseq-ids spids2 -split-sequences

    echo "${gnomon_filtering_scores_file.join('\n')}" > scores.mft
    echo "${naming_db_secondary_support.join('\n')}" > naming_db_secondary_support.mft
    echo "${search_set_secondary_support.join('\n')}" > search_set_secondary_support.mft
    echo "${best_naming_hits.join('\n')}" > best_naming_hits.mft
    
    gnomon_filter_models -models $models -nogenbank $params -quality $gnomon_quality_report -asn-cache ./asncache/  -o scored.models.asn \
            -scores-manifest scores.mft -prot-naming-hits best_naming_hits.mft  \
            -search-set-secondary-support search_set_secondary_support.mft -naming-db-secondary-support naming_db_secondary_support.mft 
    """
    stub:
    """
    touch scored.models.asn
    
    """
}
