#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow bam_strandedness {
    take:
        bam_list              // list: BAM
        sra_metadata        // path: file with sra metadata
        parameters          // Map : extra parameter and parameter update
    main:
        // Not yet used but supposed to be used by rnaseq_divide_by_strandedness 
        rnaseq_divide_by_strandedness_params = merge_params("-min-aligned 1000000 -min-unambiguous 200 -min-unambiguous-pct 2 -max-unambiguous-pct 100 -percentage-threshold 98", parameters, 'rnaseq_divide_by_strandedness')
        rnaseq_divide_by_strandedness(bam_list, sra_metadata, rnaseq_divide_by_strandedness_params)
    emit:
        strandedness = rnaseq_divide_by_strandedness.out.strandedness
        stranded_runs = rnaseq_divide_by_strandedness.out.stranded_runs
        unstranded_runs = rnaseq_divide_by_strandedness.out.unstranded_runs
        all = rnaseq_divide_by_strandedness.out.all
}


process rnaseq_divide_by_strandedness {
    label 'large_disk'
    input:
        path bam_list
        path metadata_file
        val parameters
    output:
        path "output/run.strandedness", emit: 'strandedness'
        path "output/stranded.list", emit: 'stranded_runs', optional: true
        path "output/unstranded.list", emit: 'unstranded_runs', optional: true
        path "output/*", emit: 'all' 
    script:
    """
    mkdir -p output
    mkdir -p tmp
    samtools=\$(which samtools)
    echo "${bam_list.join('\n')}" > bam_list.mft
    TMPDIR=tmp rnaseq_divide_by_strandedness -align-manifest bam_list.mft -metadata $metadata_file  $parameters  -samtools-executable \$samtools -stranded-output output/stranded.list -strandedness-output output/run.strandedness -unstranded-output output/unstranded.list
    """
    stub:
    """
    mkdir -p output
    touch output/run.strandedness
    touch output/stranded.list
    touch output/unstranded.list
    """
}
