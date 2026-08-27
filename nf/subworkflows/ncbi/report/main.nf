#!/usr/bin/env nextflow
// rnaseq long EGAPx execution
// route data to tasks

nextflow.enable.dsl=2

params.import_prefix = "../../../../nf/subworkflows/ncbi/" // redirected during testing

include { multiqc } from "./${params.import_prefix}report/multiqc/main"



workflow report_plane {
    take:
        star_logs
        busco_log_file
        mask_assm_stats
        feature_counts_xml
        feature_stats_xml
        rnaseq_long_align_report_xml
        prot_align_stats
        rnaseq_short_align_report_xml
        task_params     // task parameters for every task
    main:
        multiqc(star_logs, busco_log_file, mask_assm_stats, feature_counts_xml, feature_stats_xml, rnaseq_long_align_report_xml, prot_align_stats, rnaseq_short_align_report_xml, task_params.get('multiqc', [:]))
    emit:
        multiqc_report = multiqc.out.multiqc_report
}
