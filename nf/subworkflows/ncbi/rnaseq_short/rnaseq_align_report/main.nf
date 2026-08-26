#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow rnaseq_align_report {
    take:
        gencoll_asn
        input_metadata
        intron_counts
        run_stats
        run_list
        parameters  // Map : extra parameter and parameter update
    main:
        String report_params =  merge_params("-tracking-server NONE", parameters, 'rnaseq_align_report')
        run_align_report(gencoll_asn, input_metadata, intron_counts, run_stats, run_list, report_params)
    emit:
        align_report = run_align_report.out.align_report
        run_reports = run_align_report.out.run_reports
}


process run_align_report {
    label 'single_cpu'
    label 'small_mem'
    input:
        path gencoll_asn
        path input_metadata
        path intron_counts
        path run_stats
        val run_list
        val params
    output:
        path "rnaseq_align_report.xml", emit: "align_report"
        path "*_runs.txt", emit: "run_reports"
    script:
    """
    echo "${input_metadata.join('\n')}" > ./input_metadata.mft
    echo "${intron_counts.join('\n')}" > ./intron_counts.mft
    echo "${run_stats.join('\n')}" > ./run_stats.mft
    echo "${run_list.join('\n')}" > ./run_list.txt

    rnaseq_align_report -include-pub-and-sample-name no -gencoll-asn ${gencoll_asn} -input-manifest ./input_metadata.mft -intron-counts-manifest ./intron_counts.mft -run-stats-manifest ./run_stats.mft -run-list ./run_list.txt $params -output ./rnaseq_align_report.xml -run-report-output '@-rnaseq_runs.txt'
    # rnaseq_align_report can fail if there are no runs, so we stub the output files to avoid workflow failure
    touch rnaseq_align_report.xml
    if [ ! -s *_rnaseq_runs.txt ]; then
        touch rnaseq_runs.txt
    fi
    """
    stub:
    """
    touch rnaseq_align_report.xml
    touch stub_report_runs.txt
    """
}

