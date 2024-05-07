#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


process exec {
    label 'huge_job'
    label 'long_job'
    input:
        path seqid_list
        tuple val(sampleID), path(fasta_rna_file)
        path genome_file, stageAs: 'genome/*'
        path Star_Index
        val parameters
    output:
        path "${assembly}-${sampleID}-Aligned.out.Sorted.bam", emit: 'align'
        path "${assembly}-${sampleID}-Aligned.out.Sorted.bam.bai", emit: 'align_index'
        // path "per_run_counts.txt", emit: 'per_run_counts'
    script:
        assembly=genome_file.baseName.toString().replaceFirst(/\.(fa(sta)?|asn[bt]?)$/, "")
        if ( fasta_rna_file[0] && fasta_rna_file[1] ) {
            query_str = "${fasta_rna_file[0]},${fasta_rna_file[1]}"
        } else {
            query_str = fasta_rna_file[0]
        }
    """
    echo "${seqid_list.join('\n')}" > seqid_list.mft
    lds2_indexer -source genome
    mkdir -p out
    mkdir -p wrkarea
    #echo "${fasta_rna_file[0]}"
    #echo "${fasta_rna_file[1]}"
    echo "$query_str"
    #echo "<job query =\\\"lcl|${fasta_rna_file[0]}\\\" subject=\\\"$Star_Index\\\"></job>"
    ##echo "<job query =\\\"lcl|${fasta_rna_file[0]},${fasta_rna_file[1]}\\\" subject=\\\"$Star_Index\\\"></job>" > jobfile
    echo "<job query =\\\"lcl|${query_str}\\\" subject=\\\"$Star_Index\\\"></job>" > jobfile
    star=\$(which star-with-filter)
    samtools=\$(which samtools)
    fastq=\$(which fasterq-dump)
    star_wnode $parameters -input-jobs jobfile -genome-sequences-manifest seqid_list.mft -fastq-executable \$fastq  -samtools-executable \$samtools -star-executable \$star -output-dir . -work-area wrkarea -O out -lds2 genome/lds2.db
    """
    
    stub:
        assembly=genome_file.baseName.toString().replaceFirst(/\.(fa(sta)?|asn[bt]?)$/, "")
    """
    touch ${assembly}-${sampleID}-Aligned.out.Sorted.bam
    touch ${assembly}-${sampleID}-Aligned.out.Sorted.bam.bai
    """
}


workflow star_wnode {
    take: 
        seqid_list  // path:  list of seq ids in the index (SEQID_LIST)
        reads       // channel: FASTA file pairs generated from SRA reads, see e.g., https://www.nextflow.io/docs/latest/channel.html#fromfilepairs
        genome_file //path: genome file, fasta or ASN
        star_index  // path: index path
        parameters  // Map : extra parameter and parameter update
    main:
   
        star_params = "\" --alignSJoverhangMin 8 --outFilterMultimapNmax 200 --outFilterMismatchNmax 50 --runThreadN 16 --genomeLoad NoSharedMemory --outSAMtype SAM --outSAMattributes NH HI AS nM NM MD jM jI XS MC --outSAMprimaryFlag AllBestScore --outFilterMultimapScoreRange 50 --seedSearchStartLmax 15 --limitOutSAMoneReadBytes 1000000 --outSJtype None \""
        default_params = "-cpus-per-worker 4 -csi-threshold 512000000  -max-intron 1200000  -star-params $star_params -preserve-star-logs"
        effective_params = merge_params(default_params, parameters, 'star_wnode')
        // println("Effective parameters: $effective_params")
        exec(seqid_list, reads, genome_file, star_index, effective_params)
    emit:
        align = exec.out.align
        align_index = exec.out.align_index
}


workflow star_wnode_simplified {
    take:
        seqid_list  // path:  list of seq ids in the index (SEQID_LIST)
        reads       // list of FASTA read files, expects pairs in form SRAxxx.1, SRAxxx.2
        genome_file //path: genome file, fasta or ASN
        star_index  // path: index path
        parameters  // Map : extra parameter and parameter update
    main:
        def filePairs = Channel.of(reads).flatten().map { read ->
            def m = read =~ /\/([A-Za-z0-9]+).{1,2}$/
            [m[0][1], read]
        }.groupTuple()
        star_params = "\" --alignSJoverhangMin 8 --outFilterMultimapNmax 200 --outFilterMismatchNmax 50 --runThreadN 16 --genomeLoad NoSharedMemory --outSAMtype SAM --outSAMattributes NH HI AS nM NM MD jM jI XS MC --outSAMprimaryFlag AllBestScore --outFilterMultimapScoreRange 50 --seedSearchStartLmax 15 --limitOutSAMoneReadBytes 1000000 --outSJtype None \""
        default_params = "-cpus-per-worker 4 -csi-threshold 512000000  -max-intron 1200000  -star-params $star_params -preserve-star-logs"
        effective_params = merge_params(default_params, parameters, 'star_wnode')
        exec(seqid_list, filePairs, genome_file, star_index, effective_params)
    emit:
        align = exec.out.align
        align_index = exec.out.align_index
}
