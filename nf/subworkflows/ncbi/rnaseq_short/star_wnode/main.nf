#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


def get_effective_params(parameters, use_zcat, max_intron) {
    def effective_star_params = ""
    // Check that parameters for star_wnode doesn't contain 'star-params'
    boolean back_compat = parameters.get("star_wnode", "").contains("-star-params")
    if (!back_compat) {
        def star_params = "--alignSJoverhangMin 8 --outFilterMultimapNmax 200 --outFilterMismatchNmax 50 --runThreadN 16 --genomeLoad NoSharedMemory --outSAMtype SAM --outSAMattributes 'NH HI AS nM NM MD jM jI XS MC' --outSAMprimaryFlag AllBestScore --outFilterMultimapScoreRange 50 --seedSearchStartLmax 15 --limitOutSAMoneReadBytes 1000000 --outSJtype None"
        effective_star_params = ' -star-params "' + (use_zcat ? "--readFilesCommand zcat " : "") + merge_params(star_params, parameters, 'star-params') + '"'
    }
    def default_params = "-cpus-per-worker 4 -csi-threshold 512000000 -max-intron ${max_intron} -preserve-star-logs"
    def effective_params = merge_params(default_params, parameters, 'star_wnode') + effective_star_params
    if (!back_compat) {
        // Ad-hoc post processing - remove single quotes from effective_params
        effective_params = effective_params.replaceAll("'", "")
    }
    return effective_params
}


process run_star {
    label 'huge_job'
    label 'long_job'
    input:
        path seqid_list
        tuple val(sampleID), path(fasta_rna_file)
        val  use_zcat
        path genome_file, stageAs: 'genome/*'
        path Star_Index
        val  max_intron
        val  parameters
    output:
        path "*-Aligned.out.Sorted.bam", emit: 'align'
        path "*-Aligned.out.Sorted.bam.bai", emit: 'align_index'
        // path "per_run_counts.txt", emit: 'per_run_counts'
    script:
        def assembly=genome_file.baseName.toString().replaceFirst(/\.(fa(sta)?|asn[bt]?)$/, "")
        seqkit_cmd = ""
        if ( fasta_rna_file[0] && fasta_rna_file[1] ) {
            query_str = "${fasta_rna_file[0]},${fasta_rna_file[1]}"
            seqkit_cmd = "seqkit stats ${fasta_rna_file[0]} ${fasta_rna_file[1]}"
        } else {
            query_str = fasta_rna_file[0]
            seqkit_cmd = "seqkit stats ${fasta_rna_file[0]} "
        }
        def effective_params = get_effective_params(parameters, use_zcat, max_intron)
        // println("Effective STAR parameters: $effective_params")
    """
    echo "Assembly: ${assembly} sampleID: ${sampleID} Query: ${query_str}"
    echo "${seqid_list.join('\n')}" > seqid_list.mft
    lds2_indexer -source genome
    mkdir -p out
    mkdir -p wrkarea
    echo "<job query =\\\"lcl|${query_str}\\\" subject=\\\"$Star_Index\\\"></job>" > jobfile
    star=\$(which star-with-filter)
    samtools=\$(which samtools)
    fastq=\$(which fasterq-dump)
    ${seqkit_cmd};
    star_wnode ${effective_params} -input-jobs jobfile -genome-sequences-manifest seqid_list.mft -fastq-executable \$fastq  -samtools-executable \$samtools -star-executable \$star -output-dir . -work-area wrkarea -O out -lds2 genome/lds2.db
    """
    
    stub:
        def assembly=genome_file.baseName.toString().replaceFirst(/\.(fa(sta)?|asn[bt]?)$/, "")
        println("Assembly: ${assembly} sampleID: ${sampleID}, max_intron: ${max_intron}")
        def effective_params = get_effective_params(parameters, use_zcat, max_intron)
        println("Effective STAR parameters: $effective_params")
    """
    touch ${assembly}-${sampleID}-Aligned.out.Sorted.bam
    touch ${assembly}-${sampleID}-Aligned.out.Sorted.bam.bai
    """
}


workflow star_wnode {
    take: 
        seqid_list  // path:  list of seq ids in the index (SEQID_LIST)
        reads       // channel: FASTA file pairs generated from SRA reads, see e.g., https://www.nextflow.io/docs/latest/channel.html#fromfilepairs
        genome_file // path: genome file, fasta or ASN
        star_index  // path: index path
        max_intron  // int: maximum intron length
        parameters  // Map : extra parameter and parameter update
    main:
   
        def use_zcat_ch = reads.map { it[1][0] ==~ /.*gz$/ }
        run_star(seqid_list, reads, use_zcat_ch, genome_file, star_index, max_intron, parameters)
    emit:
        align = run_star.out.align
        align_index = run_star.out.align_index
}


workflow star_wnode_simplified {
    take:
        seqid_list  // path:  list of seq ids in the index (SEQID_LIST)
        reads       // list of FASTA read files, expects pairs in form SRAxxx.1, SRAxxx.2
        genome_file // path: genome file, fasta or ASN
        star_index  // path: index path
        max_intron  // int: maximum intron length
        parameters  // Map : extra parameter and parameter update
    main:
        def filePairs = Channel.of(reads).flatten().map { read ->
            def m = read =~ /\/([^\/]+).[12](.gz)?$/
            [m[0][1], read]
        }.groupTuple()
        def use_zcat_ch = filePairs.map { it[1][0] ==~ /.*gz$/ }
        run_star(seqid_list, filePairs, use_zcat_ch, genome_file, star_index, max_intron, parameters)
    emit:
        align = exec.out.align
        align_index = exec.out.align_index
}
