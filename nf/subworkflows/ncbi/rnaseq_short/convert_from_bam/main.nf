#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow bam2asn {
    take:
        ch_bam          // channel: BAM
        strandedness    // path: file with strandedness info for assemblies in BAM files
        genome          // path: genome in either FASTA or ASN format
        parameters      // Map: extra parameter and parameter update
    main:
        conv_param = "-filter 'pct_identity_gap >= 85' -ofmt seq-align-compressed -collapse-identical -no-scores -ifmt bam"
        convert(ch_bam, strandedness, genome, merge_params(conv_param, parameters, 'sam2asn'))
    emit:
        align = convert.out.align
        keylist = convert.out.keylist
}


process convert {
    input:
        path in_bam
        path strandedness
        path genome, stageAs: 'genome/*'
        val  conv_param
    output:
        path "${prefix}.align.asnb.gz", emit: 'align'
        path "${prefix}.align_counts.txt", emit: 'keylist'
    script:
        prefix = in_bam.name.replaceAll(/\.bam$/, '')
    """
    tmpdir=`mktemp -d`
    samtools=`which samtools`
    lds2_indexer -source genome/ -db LDS2
    # EXCEPTION_STACK_TRACE_LEVEL=Warning DEBUG_STACK_TRACE_LEVEL=Warning DIAG_POST_LEVEL=Trace
    sam2asn $conv_param -refs-local-by-default  -nogenbank -lds2 LDS2 -tmp-dir \$tmpdir -align-counts "${prefix}.align_counts.txt" -o "${prefix}.align.asnb.gz" -strandedness $strandedness -input $in_bam -samtools-path \$samtools
    """
}
