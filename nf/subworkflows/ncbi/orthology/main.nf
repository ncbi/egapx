#!/usr/bin/env nextflow
// main nextflow script for EGAPx execution
// route data to subworkflows

nextflow.enable.dsl=2

params.import_prefix = "../../../../nf/subworkflows/ncbi/" // redirected during testing

include { extract_products_from_models } from "${params.import_prefix}orthology/extract_products_from_models/main"
include { find_orthologs;} from "${params.import_prefix}orthology/find_orthologs/main"
include { fetch_ortholog_references; } from "./find_orthologs/main"
include { diamond_orthology } from "${params.import_prefix}orthology/diamond_orthology/main"
include { setup_ext_genome; setup_ext_proteins } from './../setup/main'
include { get_swiss_prot_ids as get_prot_ref_ids } from '../shared/diamond/main'

params.intermediate = false


workflow orthology_plane {
    take:
        genome_asnb
        gencoll_asn
        models
        annot_files
        ortho_files     // ref side input files
        task_params     // task parameters for every task
    main:
        // Protein alignments
        //fetch_ortholog_references()
        //fetch_ortholog_references(ortho_files)
        fetch_ortholog_references(ortho_files['genomic.fna'], ortho_files['genomic.gff'], ortho_files['protein.faa'], ortho_files['name_from.rpt'])
        def (scaffolds_ref, gencoll_ref_asn, unpacked_genome_ref, genome_ref_asn, genome_ref_asnb) = setup_ext_genome(fetch_ortholog_references.out.ref_genomic_fna, [], task_params.get('setup', [:]))
        def (unpacked_proteins_ref, proteins_ref_asn, proteins_ref_asnb) = setup_ext_proteins(fetch_ortholog_references.out.ref_protein_faa, task_params.get('setup', [:]))
        def prot_ref_ids  = get_prot_ref_ids(proteins_ref_asnb)
        //orthology plane
        extract_products_from_models(annot_files, task_params.get('extract_products_from_models', [:]))
        // reference side. 1) gencoll, annotation , genome and protein sequence --> asn cashe or LDS demon
        
        diamond_orthology(extract_products_from_models.out.prot_ids, prot_ref_ids , models, proteins_ref_asnb , task_params.get('diamond_orthology', [:]))

        // input side 1) gencoll asn from setup, genome_asn from setup, protein from gnomon_wnode.out, annotation it is from annotbuilder annot_files or accepts_asn
        find_orthologs( gencoll_asn,  gencoll_ref_asn,  annot_files, fetch_ortholog_references.out.annot_file, diamond_orthology.out.alignments, [],
                        models, proteins_ref_asnb , genome_asnb, genome_ref_asnb ,  task_params.get('find_orthologs', [:]))

    emit:
        orthologs = find_orthologs.out.orthologs
        name_from_ortholog = fetch_ortholog_references.out.name_from_ortholog
}

