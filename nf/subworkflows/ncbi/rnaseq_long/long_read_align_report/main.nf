#!/usr/bin/env nextflow
nextflow.enable.dsl=2
include { merge_params } from '../../utilities'

workflow long_read_align_report {
    take:
        gencoll_asn
        minimap_stats
        sra_metadata
        filtered_stats
        run_id_map
        parameters  // Map : extra parameter and parameter update
    main:
        String long_read_align_report_params = merge_params('-tracking-server NONE -tracking-user Username -tracking-password Locator', parameters, 'long_read_align_report')
        run_long_read_align_report(gencoll_asn, minimap_stats, sra_metadata, filtered_stats, run_id_map, long_read_align_report_params)
    emit:
        align_report = run_long_read_align_report.out.align_report
}

process run_long_read_align_report {
    label 'single_cpu'
    label 'small_mem'
    input:
        path gencoll_asn
        path minimap_stats
        path sra_metadata
        path filtered_stats
        path run_id_map
        val parameters
    output:
        path "output/long_read_report.xml", emit: 'align_report'
    script:
    """
    mkdir -p output
    # adding logic to filter out relevant rows from sra_metadata.dat, and using that for the manifest
    grep "^#" ${sra_metadata} > sra_metadata_filtered.dat
    grep "^#" -v ${minimap_stats} | cut -f 2 | fgrep -f - ${sra_metadata} >> sra_metadata_filtered.dat
    echo "sra_metadata_filtered.dat" > sra_metadata.mft
    echo "${minimap_stats.join('\n')}" > minimap_stats.mft
    echo "${filtered_stats.join('\n')}" > filtered_stats.mft
    # long_read_align_report does not need -nogenbank parameter, nor does it support it.
    long_read_align_report -filtered-stats-manifest filtered_stats.mft  $parameters \
        -gencoll-asn $gencoll_asn  -minimap-stats-manifest minimap_stats.mft  -output output/long_read_report.xml  -sra-metadata-manifest sra_metadata.mft
    # Remap renamed run IDs back to original names in the output XML
    if [ -n "${run_id_map}" ] && [ -s "${run_id_map}" ]; then
        while IFS=\$'\\t' read -r renamed_id original_id; do
            sed -i "s/\${renamed_id}/\${original_id}/g" output/long_read_report.xml
        done < "${run_id_map}"
    fi
    """
    stub:
    """
    mkdir -p output
    touch output/long_read_report.xml
    """
}
