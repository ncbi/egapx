#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'

workflow locus_link {
    take:
        best_refseq_prot_hit
        orthologs
        annotation 
        gencoll_asn
        gnomon_lds2_source
        best_prot_hit
        track_loci
        comparisons
        curr_prev_compare
        gnomon_biotypes
        lxr_data
        proteins_asn
        name_from_ortholog
        parameters  // Map : extra parameter and parameter update
    main:
        default_params = ""
        effective_params = merge_params(default_params, parameters, 'locus_type')
        run_locus_link(best_refseq_prot_hit, orthologs, annotation, 
                 gencoll_asn, gnomon_lds2_source, best_prot_hit, track_loci, comparisons,  curr_prev_compare,  
                 gnomon_biotypes, lxr_data, proteins_asn, name_from_ortholog,  default_params)
    emit:
        best_gnomon_prot_hit = run_locus_link.out.best_gnomon_prot_hit
        best_refseq_prot_hit = run_locus_link.out.best_refseq_prot_hit
        locustypes = run_locus_link.out.locustypes
        locus = run_locus_link.out.locus
        stats = run_locus_link.out.stats
        all = run_locus_link.out.all
}



process run_locus_link {
    input:
        path best_refseq_prot_hit
        path orthologs
        path annotation 
        path gencoll_asn
        path gnomon_lds2_source, stageAs: 'genome/*'
        path best_prot_hit
        path track_loci
        path comparisons
        path curr_prev_compare
        path gnomon_biotypes
        path lxr_data, stageAs: 'lxr_tracking_data.txt'
        path proteins_asn
        path name_from_ortholog
        val parameters
    output:
        path ('output/best_gnomon_prot_hit.tsv'), emit: 'best_gnomon_prot_hit'
        path ('output/best_refseq_prot_hit.tsv'), emit: 'best_refseq_prot_hit'
        path ('output/locustypes.tsv'), emit: 'locustypes'
        path ('output/locus.lnk'), emit: 'locus'
        path ('output/stats.xml'), emit: 'stats'
        path ('output/*'), emit: 'all'
    script:
    """
    mkdir -p output
    mkdir -p ./asncache/
    prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i $proteins_asn  -oseq-ids /dev/null -split-sequences
    
    lds2_indexer -source genome/ -db LDS2 
    echo "${best_prot_hit.join('\n')}" > best_prot_hit.mft
    extract_prot_names -alns best_prot_hit.mft  -nogenbank -o output/best_gnomon_prot_hit.tsv -asn-cache ./asncache/ -lds2 LDS2
    echo "${best_refseq_prot_hit.join('\n')}" > best_refseq_prot_hit.mft
    extract_prot_names -alns best_refseq_prot_hit.mft -nogenbank -o output/best_refseq_prot_hit.tsv -asn-cache ./asncache/ -lds2 LDS2 
    echo "${annotation.join('\n')}" > annotation.mft
    echo "${curr_prev_compare.join('\n')}" > curr_prev_compare.mft
    echo  "${comparisons.join('\n')}"  > comparisons.mft
    str=""
    str="\$str -orthologs $orthologs"
    str="\$str -lxr $lxr_data"
    str="\$str -locus_track $track_loci"
    str="\$str -name_from_ortholog_rpt $name_from_ortholog"

    locus_type  -asn-cache ./asncache/ -lds2 ./LDS2 -nogenbank -no_acc_reserve  -annots annotation.mft -gc $gencoll_asn -gnomon_biotype $gnomon_biotypes -o_stats output/stats.xml -o_locustypes output/locustypes.tsv -o_locus_lnk output/locus.lnk  -annotcmp comparisons.mft  -annotcmp_pb curr_prev_compare.mft \$str
    """
    stub:
    """
    mkdir -p output
    touch output/best_gnomon_prot_hit.tsv
    touch output/best_refseq_prot_hit.tsv
    touch output/locustypes.tsv
    touch output/locus.lnk
    touch output/stats.xml
    """
}

