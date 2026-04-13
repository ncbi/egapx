#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow bam2asn {
    take:
        ch_bam          // channel: BAM
        strandedness    // path: file with strandedness info for assemblies in BAM files
        genome          // path: genome in either FASTA or ASN format
        sra_metadata    // path: sra metadata tsv file
        parameters      // Map: extra parameter and parameter update
    main:
        conv_param = "-filter 'pct_identity_gap >= 85' -ofmt seq-align-compressed -collapse-identical -no-scores -ifmt bam"
        convert(ch_bam, strandedness, genome, sra_metadata, merge_params(conv_param, parameters, 'sam2asn'))
    emit:
        align = convert.out.align
        keylist = convert.out.keylist
}


process convert {
    label 'long_job'
    label 'large_disk'
    label 'multi_node'
    label 'small_mem'
    input:
        path in_bam
        path strandedness
        path genome, stageAs: 'genome/*'
        path sra_metadata
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
    mkdir -p tmp/sam
    if [[ -n \${TMPDIR-} ]]; then
        mkdir -p \${TMPDIR} || true
    fi
    if [[ -n \${TEMP-} ]]; then
        mkdir -p \${TEMP} || true
    fi

    lds2_indexer -source genome/ -db tmp/LDS2
    sam2asn $conv_param -pseudo-run-accessions $sra_metadata -refs-local-by-default  -nogenbank -lds2 tmp/LDS2 -tmp-dir tmp/sam -align-counts "${prefix}.align_counts.txt" -o "${prefix}.align.asnb.gz" -strandedness $strandedness -input $in_bam -samtools-path \$samtools
    rm -rf
    """

    stub:
        prefix = in_bam.name.replaceAll(/\.bam$/, '')
    """
    touch ${prefix}.align_counts.txt
    touch ${prefix}.align.asnb.gz
    """
}
