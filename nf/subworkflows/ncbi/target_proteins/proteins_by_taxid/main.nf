#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow proteins_by_taxid {
    take:
        genome_asn     //path: genome file
        proteins_asn   //path: protein file
        tab_file       //path: alignment asn file from paf2asn
        parameters     // Map : extra parameter and parameter update
    main:
        String param = merge_params("", parameters, "products_by_taxid")
        run_proteins_by_taxid(genome_asn, proteins_asn, tab_file, param)
    emit:
        taxid_list = run_proteins_by_taxid.out.taxid_list
        id_files = run_proteins_by_taxid.out.id_files
        mft_files = run_proteins_by_taxid.out.mft_files
}


process run_proteins_by_taxid{
    label 'single_cpu'
    label 'small_mem'
    input:
        path genome, stageAs: 'indexed/genome.asnt'
        path proteins,  stageAs: 'indexed/proteins.asnt'
        path tab_file
        val parameters
    output:
        path "output/taxid.list", emit: 'taxid_list'
        path "output/*.ids", emit: 'id_files'
        path "output/*.mft", emit: 'mft_files'
    script:
    """
    mkdir -p output
    mkdir -p tmp/asncache
    #lds2_indexer -source indexed -db tmp/lds_index
    products_by_taxid  $parameters -products $tab_file -products-output output/@.ids -products-output-manifest output/@.mft -taxids-output output/taxid.list
    rm -rf tmp
    """
    stub:
    """
    mkdir -p output
    touch output/taxid.list
    touch output/high_ident_xps.ids
    touch output/high_ident_xps.mft
    touch output/low_ident_xps.ids
    touch output/low_ident_xps.mft
    """
}
