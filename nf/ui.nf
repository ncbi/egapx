#!/usr/bin/env nextflow
// main nextflow script for EGAPx ui execution
// prepare data channels and call main subworkflow

nextflow.enable.dsl=2

include { egapx } from './subworkflows/ncbi/main'


process export {
    publishDir "${params.output}", mode: 'copy', saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    input:
        path genomic_gff
        path genomic_gtf
        path genomic_fasta
        path transcript_fasta
        path cds_fasta
        path proteins_fasta
        path annot_builder_output, stageAs: 'annot_builder_output/*'
        path validated, stageAs: 'validated/*'
        path stats, stageAs: 'stats/*'
        path annotated_genome_asn
        path annotation_data_comment
        path star_bam, stageAs: 'STAR/*'
        path busco_results, stageAs: 'busco/*'
        path gnomon_biotype_contam_rpt, stageAs: 'GNOMON/*'
        path user_gnomon_report, stageAs: 'GNOMON/*'
        path user_gnomon_quality_report, stageAs: 'GNOMON/*'
        path gnomon_summaries, stageAs: 'GNOMON/*'
        // path locus
    output:
        path "*", includeInputs: true
    script:
    """
    echo "export script"
    """
    stub:
    """
    echo "export stub"
    """
}


workflow {
    // Parse input parameters
    def input_params = params.get('input', [:])
    def task_params = params.get('tasks', [:])
    egapx(input_params, task_params)
    export(egapx.out.out_ggff, 
           egapx.out.out_ggtf,
           egapx.out.out_gfa,
           egapx.out.out_rna_fa,
           egapx.out.out_cds_fa,
           egapx.out.out_prot_fa,
           egapx.out.annot_builder_output,
           egapx.out.validated,
           egapx.out.stats,
           egapx.out.annotated_genome_asn,
           egapx.out.annotation_data_comment,
           egapx.out.out_star_bam,
           egapx.out.busco_results,
           egapx.out.gnomon_biotype_contam_rpt,
           egapx.out.user_gnomon_report,
           egapx.out.user_gnomon_quality_report,
           egapx.out.gnomon_summaries)
}
