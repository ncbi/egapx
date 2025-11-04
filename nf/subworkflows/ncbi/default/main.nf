#!/usr/bin/env nextflow
// gnomon plane workflow
// route data to tasks

nextflow.enable.dsl=2

params.import_prefix = "../../../../nf/subworkflows/ncbi/" // redirected during testing


include { get_hmm_params} from "./${params.import_prefix}default/get_hmm_params/main"
include { convert_summary_files} from "${params.import_prefix}default/convert_summary_files/main"
include { convert_annotations } from "${params.import_prefix}default/convert_annotations/main"
include { align_sort_sa } from "${params.import_prefix}default/align_sort_sa/main"

//include { locus_track } from './locus_track/main'
//include { locus_link } from './locus_link/main'


params.intermediate = false

workflow default_plane {
    take:
        genome_asn
        proteins_asn
        alignments // list of all relevent input alignments
        asn_file
        report
        quality_report
        locustypes
        locus_tag_prefix

        tax_id          // NCBI tax id of the closest taxon to the genome

        task_params     // task parameters for every task
    main:
        // GNOMON
        get_hmm_params(tax_id, task_params.get('get_hmm_params', [:]))
        // Seed Protein-Model Hits
        convert_summary_files(report, quality_report, locustypes, locus_tag_prefix, task_params.get('convert_summary_files', [:]))
        convert_annotations(asn_file , task_params.get('convert_annotations', [:]))
        sorted_asn_file = align_sort_sa(genome_asn, proteins_asn, alignments, task_params.get('align_sort_sa', [:]))  

    emit:
        hmm_params = get_hmm_params.out.outputs
        gnomon_quality_report = quality_report
        gnomon_report = report
        out_ggff = convert_annotations.out.genomic_gff
        out_ggtf = convert_annotations.out.genomic_gtf
        out_gfa  = convert_annotations.out.genomic_fasta
        out_rna_fa = convert_annotations.out.transcripts_fasta
        out_cds_fa = convert_annotations.out.cds_fasta
        out_prot_fa = convert_annotations.out.proteins_fasta
        sorted_asn_file  = sorted_asn_file
}