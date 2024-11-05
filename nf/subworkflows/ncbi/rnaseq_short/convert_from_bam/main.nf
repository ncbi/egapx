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
    label 'long_job'
    label 'large_disk'
    input:
        path in_bam
        path strandedness
        path genome, stageAs: 'genome/*'
        val  conv_param
    output:
        path "${prefix}.align.asnb.gz", emit: 'align', optional: true
        path "${prefix}.align_counts.txt", emit: 'keylist', optional: true
    script:
        prefix = in_bam.name.replaceAll(/\.bam$/, '')
        min_file_size = 100000
    """
    samtools=`which samtools`
    if [ `stat -L -c%s $in_bam` -lt $min_file_size ] && [ `\$samtools view -c $in_bam` -eq 0 ]; then
        exit 0
    fi
    tmpdir=`mktemp -d --tmpdir=.`
    lds2_indexer -source genome/ -db LDS2
    # EXCEPTION_STACK_TRACE_LEVEL=Warning DEBUG_STACK_TRACE_LEVEL=Warning DIAG_POST_LEVEL=Trace
    sam2asn $conv_param -refs-local-by-default  -nogenbank -lds2 LDS2 -tmp-dir \$tmpdir -align-counts "${prefix}.align_counts.txt" -o "${prefix}.align.asnb.gz" -strandedness $strandedness -input $in_bam -samtools-path \$samtools
    """

    stub:
        prefix = in_bam.name.replaceAll(/\.bam$/, '')
    """
    touch ${prefix}.align_counts.txt
    touch ${prefix}.align.asnb.gz
    """
}
