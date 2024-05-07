#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow rnaseq_register_stats {
    take:
        gencoll_asn
        sra_metadata
        collapsed_aligns
        per_run_counts
        stranded_runs
        unstranded_runs
        parameters  // Map : extra parameter and parameter update
    main:
        String stats_params =  merge_params(" -TS NONE -TU NONE -TP NONE ", parameters, 'rnaseq_register_stats')
        register_stats(gencoll_asn, sra_metadata, collapsed_aligns, per_run_counts, stranded_runs, unstranded_runs, stats_params)
    emit:
        outputs = register_stats.out.outputs
}


process register_stats {
    input:
        path gencoll_asn
        path sra_metadata
        path collapsed_aligns
        path per_run_counts
        path stranded_runs
        path unstranded_runs
        val params
    output:
        path "run_stats.tsv", emit: "outputs"
    script:
    """
    ##echo "${gencoll_asn.join('\n')}" > ./gencoll_asn.mft
    echo "${sra_metadata.join('\n')}" > ./sra_metadata.mft
    echo "${collapsed_aligns.join('\n')}" > ./collapsed_aligns.mft
    echo "${per_run_counts.join('\n')}" > ./per_run_counts.mft
    echo "${stranded_runs.join('\n')}" > ./stranded_runs.mft
    echo "${unstranded_runs.join('\n')}" > ./unstranded_runs.mft

    rnaseq_register_stats    \
      -gencoll-asn ${gencoll_asn}           \
      -sra-metadata-manifest ./sra_metadata.mft         \
      -collapsed-aligns-manifest ./collapsed_aligns.mft \
      -per-run-counts-manifest ./per_run_counts.mft     \
      -stranded-runs-manifest ./stranded_runs.mft       \
      -unstranded-runs-manifest ./unstranded_runs.mft   \
      $params  -o ./run_stats.tsv

    ls -l 
    """

    stub:
    """
    touch run_stats.tsv
    """
}

