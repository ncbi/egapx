#!/usr/bin/env nextflow
// gnomon plane workflow
// route data to tasks

nextflow.enable.dsl=2

params.import_prefix = "../../../../nf/subworkflows/ncbi/" // redirected during testing

include { gnomon_biotype} from "${params.import_prefix}annot_proc/gnomon_biotype/main"
include { orthology_plane } from "../orthology/main"
include { locus_track } from "${params.import_prefix}annot_proc/locus_track/main"
include { locus_link } from "${params.import_prefix}annot_proc/locus_link/main"
include { annot_builder } from "${params.import_prefix}annot_proc/annot_builder/main"
include { final_asn_markup } from "${params.import_prefix}annot_proc/final_asn_markup/main"
include { generate_fasta_from_annots } from "${params.import_prefix}annot_proc/fasta_from_annots/main"
include { annotwriter } from "${params.import_prefix}annot_proc/annotwriter/main"


params.intermediate = false


// curr_tax_id|is_curated_genome|symbol_format_class|extrn_naming_authority|use_in_lxr_client|<unused by tax_attrs>|orth_source_tax_id|schema_ver
// tax_id|0|symbol_format_class|1|0|y|ortho_tax_id(zero for now?)|2
process print_fake_lxr_data {
    input:
        val tax_id
        val symbol_format_class
        val ortho_tax_id
    output:
        path('fake_lxr.tsv'), emit: 'lxr_data'
    script:
    """
        echo '_TAX_ATTRS_' > 'fake_lxr.tsv'
        echo '$tax_id|0|$symbol_format_class|1|0|y|$ortho_tax_id|2' >> 'fake_lxr.tsv'
    """
    stub:
    """
        echo '_TAX_ATTRS_' > 'fake_lxr.tsv'
        echo '9606|0|allupper|1|0|y|9606|2' >> 'fake_lxr.tsv'
    """
}

workflow annot_proc_plane {
    take:
        annotation_name_prefix
        gnomon_models
        gencoll_asn
        genome_asn 
        genome_asnb
        scaffolds

        //orthologs
        //name_from_ortholog 

        tax_id          // NCBI tax id of the closest taxon to the genome
        symbol_format_class // string for how to format gene names 
        name_cleanup_rules_file // file with name correction rules
        ortho_files      /// ortho reference input files
        alignments
        best_naming_hits
        swiss_prot_asn
        prot_denylist
        task_params     // task parameters for every task
    main:
        // Post GNOMON
        // might come its own plane     
        gnomon_biotype(gnomon_models,/*splices_file  -- constant*/ [],  prot_denylist, gencoll_asn, swiss_prot_asn, [], alignments, name_cleanup_rules_file, task_params.get('gnomon_biotype', [:]))

        annot_builder(gencoll_asn, gnomon_models, genome_asn, task_params.get('annot_builder', [:]))
        def accept_ftable_file = annot_builder.out.accept_ftable_annot
        def annot_files = annot_builder.out.annot_files

        lxr_data = print_fake_lxr_data(tax_id, symbol_format_class, ortho_files.get('taxid',0)).lxr_data

        orthology_plane(genome_asnb, gencoll_asn, gnomon_models, annot_files, ortho_files, task_params)
        def orthologs = orthology_plane.out.orthologs
        def name_from_ortholog = orthology_plane.out.name_from_ortholog


        ///def lxr_data = []
        locus_track(accept_ftable_file, gencoll_asn, lxr_data, task_params.get('locus_track', [:]) )
        
        locus_link(/*best_refseq_prot_hit  -- best protein hits from refseq plane*/ [], orthologs, accept_ftable_file, 
                gencoll_asn, gnomon_models, best_naming_hits , locus_track.out.track_rpt, /*comparisons*/ [],  /*curr_prev_compare*/ [], 
                gnomon_biotype.out.biotypes, lxr_data, name_cleanup_rules_file, swiss_prot_asn, name_from_ortholog,  task_params.get('locus_link', [:]))

        // Replace spaces with underscores in assembly name
        def assembly_name = annotation_name_prefix.replaceAll(' ', '_')
   
        final_asn_markup(assembly_name, gencoll_asn, genome_asn, scaffolds, /*chromosomes ASN*/ [], annot_builder.out.accept_asn.collect(), locus_link.out.locus, locus_link.out.locustypes, task_params.get('final_asn_markup', [:]) )

        generate_fasta_from_annots(gencoll_asn, final_asn_markup.out.to_convert, /*chromosomes ASN*/ [], task_params.get('generate_fasta_from_annots', [:]))

        annotwriter(accept_ftable_file, task_params.get('annotwriter', [:]))

    emit:
        locus = locus_link.out.locus
        locustypes = locus_link.out.locustypes
        accept_annot_file = accept_ftable_file
        final_asn_out = final_asn_markup.out.outputs
        to_convert = final_asn_markup.out.to_convert
        validated = final_asn_markup.out.validated
        stats = final_asn_markup.out.stats
        annotated_genome_asn = final_asn_markup.out.annotated_genome_asn
        annotation_data_comment = final_asn_markup.out.annotation_data_comment
        gff_annotated_file = annotwriter.out.annoted_file
        annot_proteins = generate_fasta_from_annots.out.proteins
        gnomon_biotype_contam_rpt = gnomon_biotype.out.contam_rpt
}
