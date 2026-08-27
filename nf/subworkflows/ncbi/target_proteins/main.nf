#!/usr/bin/env nextflow
// main nextflow script for EGAPx execution
// route data to subworkflows

nextflow.enable.dsl=2

params.import_prefix = "../../../../nf/subworkflows/ncbi/" // redirected during testing

include { miniprot } from "./${params.import_prefix}target_proteins/miniprot/main"
include { paf2asn } from "./${params.import_prefix}target_proteins/paf2asn/main"
include { tblastn_align } from "./${params.import_prefix}target_proteins/tblastn_align/main"
include { prosplign_prepare } from "./${params.import_prefix}target_proteins/prosplign_prepare/main"
include { prosplign_wnode } from "./${params.import_prefix}target_proteins/prosplign_wnode/main"
include { best_aligned_prot } from "./${params.import_prefix}target_proteins/best_aligned_prot/main"
include { align_filter_sa } from "./${params.import_prefix}target_proteins/align_filter_sa/main"
include { align_filter_sa as align_filter_sa_remove_low_quality } from "./${params.import_prefix}target_proteins/align_filter_sa/main"
include { align_sort_sa} from "./${params.import_prefix}target_proteins/../default/align_sort_sa/main"
include { align_as_tabular_sa } from "./${params.import_prefix}target_proteins/align_as_tabular_sa/main"
include { proteins_by_taxid } from "./${params.import_prefix}target_proteins/proteins_by_taxid/main"
include { prot_align_stats } from "./${params.import_prefix}target_proteins/prot_align_stats/main"

include { checkpoint_save; checkpoint_load_channels; CHECKPOINT_LOAD_JSON} from "./../../../../nf/lib/checkpoint_both"

params.intermediate = false

workflow target_proteins_plane_cached {
    take:
        unpacked_genome_fasta
        genome_asn
        genome_blastdb
        gencoll_asn
        unpacked_proteins_fasta
        proteins_asn
        aligner_name
        max_intron
        assembly_taxid_list
        proteins_accessions_list
        proteins_filter_taxons
        annotation_name_prefix
        task_params     // task parameters for every task.
    main:
        def checkpoint_dir = task_params.get('checkpoints_dir', [])
        def checkpoint_name = 'target_proteins_plane'
        def checkpoint_enabled = task_params.get('checkpoints_save', false)
        def checkpoint_args = [name: checkpoint_name, dir: checkpoint_dir, enabled: checkpoint_enabled]
        def ck_file = file("${checkpoint_dir}/${checkpoint_name}.json")
        def out_obj=[:]
        if (ck_file.exists()) {
            CHECKPOINT_LOAD_JSON(checkpoint_name, checkpoint_dir)
            //def loaded = checkpoint_load_channels(CHECKPOINT_LOAD_JSON.out.json_file)
            def loaded = checkpoint_load_channels(checkpoint_args)
            out_obj = loaded
        } else {
            target_proteins_plane(unpacked_genome_fasta, genome_asn, genome_blastdb, 
                                gencoll_asn, unpacked_proteins_fasta, proteins_asn, 
                                aligner_name, max_intron, 
                                assembly_taxid_list, proteins_accessions_list, proteins_filter_taxons, 
                                annotation_name_prefix, task_params)
            out_obj = target_proteins_plane.out
            checkpoint_save([target_proteins_plane.out], checkpoint_args)
        }
    emit:
        protein_alignments = out_obj.protein_alignments
        filtered_protein_alignments = out_obj.filtered_protein_alignments
        prot_align_stats = out_obj.prot_align_stats
}
workflow target_proteins_plane {
    take:
        unpacked_genome_fasta
        genome_asn
        genome_blastdb
        gencoll_asn
        unpacked_proteins_fasta
        proteins_asn
        aligner_name
        max_intron
        assembly_taxid_list
        proteins_accessions_list
        proteins_filter_taxons
        annotation_name_prefix
        task_params     // task parameters for every task
    main:
        // Protein alignments
       
        alignments_to_use = []
        if (aligner_name == "miniprot" ) {
            miniprot(unpacked_genome_fasta, unpacked_proteins_fasta, max_intron, task_params.get('miniprot', [:]))
            def miniprot_file = miniprot.out.miniprot_file
            paf2asn(genome_asn, proteins_asn, miniprot_file, task_params.get('paf2asn', [:]))
            alignments_to_use = paf2asn.out.asn_file.collect(sort: true)
        } 
        else if (aligner_name == "prosplign") {
            tblastn_align(genome_asn, proteins_asn, genome_blastdb, task_params.get('tblastn_align', [:]))
            tblastn_aligns = tblastn_align.out.blast_asn.collect()
            prosplign_prepare(genome_asn, proteins_asn, tblastn_aligns, gencoll_asn, max_intron, task_params.get('prosplign_prepare', [:]))
            prosplign_wnode(genome_asn, proteins_asn, prosplign_prepare.out.compartments_asn, max_intron, task_params.get('prosplign_wnode', [:]))
            alignments_to_use = prosplign_wnode.out.prosplign_asn.collect(sort: true)
        } else {
            error "Not implemented"
        }

        best_aligned_prot(genome_asn, proteins_asn, alignments_to_use, gencoll_asn, task_params.get('best_aligned_prot', [:]))
        align_filter_sa(genome_asn, proteins_asn, best_aligned_prot.out.asn_file, [],  task_params.get('align_filter_sa', [:]))
        align_sort_sa(genome_asn, proteins_asn,align_filter_sa.out.filtered_file, task_params.get('align_sort_sa', [:]))
        align_as_tabular_sa(genome_asn, proteins_asn, align_sort_sa.out.sorted_asn_file, task_params.get('align_as_tabular_sa', [:]))
        proteins_by_taxid(genome_asn, proteins_asn, align_as_tabular_sa.out.align_tab, task_params.get('proteins_by_taxid', [:]))
        align_filter_sa_remove_low_quality(genome_asn, proteins_asn, align_sort_sa.out.sorted_asn_file, proteins_by_taxid.out.id_files, task_params.get('align_filter_sa_remove_low_quality', [:]))
        prot_align_stats(align_as_tabular_sa.out.align_tab, assembly_taxid_list, proteins_accessions_list, proteins_filter_taxons, annotation_name_prefix, task_params.get('prot_align_stats', [:])) 
    emit:
        protein_alignments = align_sort_sa.out
        filtered_protein_alignments = align_filter_sa_remove_low_quality.out.filtered_file
        prot_align_stats = prot_align_stats.out.prot_align_stats
}
