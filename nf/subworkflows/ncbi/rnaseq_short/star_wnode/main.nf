#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'

params.verbose = false


def get_effective_params(parameters, max_intron, cpus) {
    def effective_star_params = ""
    // Check that parameters for star_wnode doesn't contain 'star-params'
    boolean back_compat = parameters.get("star_wnode", "").contains("-star-params")
    if (!back_compat) {
        def star_params = "--alignSJoverhangMin 8 --outFilterMultimapNmax 200 --outFilterMismatchNmax 50 --runThreadN ${cpus} --genomeLoad NoSharedMemory --outSAMtype SAM --outSAMattributes 'NH HI AS nM NM MD jM jI XS MC' --outSAMprimaryFlag AllBestScore --outFilterMultimapScoreRange 50 --seedSearchStartLmax 15 --limitOutSAMoneReadBytes 1000000 --outSJtype None"
        effective_star_params = ' -star-params "--readFilesCommand zstdcat ' + merge_params(star_params, parameters, 'star-params') + '"'
    }
    def default_params = "-cpus-per-worker ${cpus} -csi-threshold 0 -max-intron ${max_intron} -preserve-star-logs"
    def effective_params = merge_params(default_params, parameters, 'star_wnode') + effective_star_params
    // Force CSI threshold to 0 to always use CSI
    effective_params = effective_params.replaceAll(/-csi-threshold [0-9]+/, "-csi-threshold 0")
    if (!back_compat) {
        // Ad-hoc post processing - remove single quotes from effective_params
        effective_params = effective_params.replaceAll("'", "")
    }
    return effective_params
}


process run_star {
    label 'multi_node'
    label 'med_mem'
    label 'long_job'
    input:
        path seqid_list
        tuple val(sampleID), path(fasta_rna_file)
        path genome_file, stageAs: 'genome/*'
        path Star_Index
        val  max_intron
        val  parameters
    output:
        path "*-Aligned.out.Sorted.bam", emit: 'align'
        path "*-Aligned.out.Sorted.bam.csi", emit: 'align_index'
        // path "per_run_counts.txt", emit: 'per_run_counts'
    script:
        def assembly=genome_file.baseName.toString().replaceFirst(/\.(fa(sta)?|asn[bt]?)$/, "")
        seqkit_cmd = ""
        if ( fasta_rna_file[0] && fasta_rna_file[1] ) {
            query_str = "${fasta_rna_file[0]},${fasta_rna_file[1]}"
            seqkit_cmd = "echo ${fasta_rna_file[0]}; zstdcat ${fasta_rna_file[0]} | seqkit stats\necho ${fasta_rna_file[1]}; zstdcat ${fasta_rna_file[1]} | seqkit stats"
        } else {
            query_str = fasta_rna_file[0]
            seqkit_cmd = "echo ${fasta_rna_file[0]}; zstdcat ${fasta_rna_file[0]} | seqkit stats"
        }

        // For some executors (e.g. SGE) the task.cpus is not set correctly, they allocate correct number of threads through clusterOptions.
        // We use ext.cpus to pass the number of cpus here and use clusterOptions to allocate large enough instance
        def effective_params = get_effective_params(parameters, max_intron, task.ext.threads)
        if (params.verbose) {
            println("Effective STAR parameters: $effective_params")
        }
    """
    echo "Assembly: ${assembly} sampleID: ${sampleID} Query: ${query_str}"
    echo "${seqid_list.join('\n')}" > seqid_list.mft
    lds2_indexer -source genome
    mkdir -p out
    mkdir -p wrkarea
    if [[ -n \${TMPDIR-} ]]; then
        mkdir -p \${TMPDIR} || true
    fi
    if [[ -n \${TEMP-} ]]; then
        mkdir -p \${TEMP} || true
    fi
    echo "<job query =\\\"lcl|${query_str}\\\" subject=\\\"$Star_Index\\\"></job>" > jobfile
    star=\$(which star-with-filter)
    samtools=\$(which samtools)
    fastq=\$(which fasterq-dump)
    seqkit=\$(which seqkit)
    ${seqkit_cmd}

    star_wnode ${effective_params} -input-jobs jobfile -genome-sequences-manifest seqid_list.mft -seqkit-executable \$seqkit -fastq-executable \$fastq  -samtools-executable \$samtools -star-executable \$star -output-dir . -work-area wrkarea -O out -lds2 genome/lds2.db

    # re-header BAM, dropping @PG and @CO lines containing non-deterministic elements.
    for bam in *.bam; do
        hdr="\$bam".hdr
        samtools view -H "\$bam" | grep -vE '^@PG|^@CO' > "\$hdr"
        samtools reheader --no-PG "\$hdr" "\$bam" > "\$bam".new  # NB: couldn't get --in-place to work
        rm -f "\$hdr"
        mv "\$bam".new "\$bam"

        # recreate index
        [[ -f "\${bam}.bai" ]] && samtools index -@ 4 -b "\$bam"
        [[ -f "\${bam}.csi" ]] && samtools index -@ 4 -c "\$bam"
    done
    """
    
    stub:
        def assembly=genome_file.baseName.toString().replaceFirst(/\.(fa(sta)?|asn[bt]?)$/, "")

        println("Assembly: ${assembly} sampleID: ${sampleID}, max_intron: ${max_intron}")
        def effective_params = get_effective_params(parameters, max_intron, task.ext.threads)
        println("Effective STAR parameters: $effective_params")
    """
    # NB: see GP-40504
    echo ${task.index} > ${assembly}-${sampleID}-Aligned.out.Sorted.bam
    echo ${task.index} > ${assembly}-${sampleID}-Aligned.out.Sorted.bam.csi
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
        run_star(seqid_list, reads, genome_file, star_index, max_intron, parameters)
    emit:
        align = run_star.out.align
        align_index = run_star.out.align_index
}
