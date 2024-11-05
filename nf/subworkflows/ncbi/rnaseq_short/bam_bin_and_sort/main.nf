#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow bam_bin_and_sort {
    take:
        ch_bam          // channel: BAM
        ch_index        // channel: BAM index
        genome          // path: genome fasta
        organelle       // list of path: organelles
        parameters      // Map : extra parameter and parameter update
    main:
        def assembly_sizes = calc_assembly_sizes(ch_bam.collect())
        def bam_bin_params = " -file-pattern 'bin#.bam' "
        out = bam_bin(ch_bam, ch_index, genome, organelle, assembly_sizes, merge_params(bam_bin_params, parameters, 'bam_bin'))
        bins = out.collect()
        merge_args = merge_prepare(bins)
        // HACK: relies on the same filenames for files in bins as
        // they were at the merge_prepare step.
        merge(merge_args.flatten(), bins)
    emit:
        sorted = merge.out
}


process calc_assembly_sizes {
    label 'large_disk'
    input:
        path bam_files
    output:
        path "assembly_sizes.hash"
    script:
    """
    #!/usr/bin/env bash

    declare -A assembly_sizes
    for s in ${bam_files}; do
        chunk_size=`wc -c \${s} | awk '{print \$1}'`
        # Select assembly name
        s=\$(basename \${s})
        regex="^(.+)-([^-]+)-Aligned[.]out[.]Sorted[.]bam\$"
        if [[ \$s =~ \$regex ]]; then
            assembly="\${BASH_REMATCH[1]}"
        else
            echo "Malformed BAM name, \${s}"
            exit 1
        fi
        # Update assembly -> total_size map
        echo "Assembly: \${assembly}, chunk size: \${chunk_size}"
        if [[ -z "\${assembly_sizes[\${assembly}]}" ]]; then
            assembly_sizes[\${assembly}]=\${chunk_size}
        else
            assembly_sizes[\${assembly}]=\$((assembly_sizes[\${assembly}] + chunk_size))
        fi
    done
    declare -p assembly_sizes > assembly_sizes.hash
    """

    stub:
    """
    touch assembly_sizes.hash
    """
}


process bam_bin {
    input:
        path sorted_bam
        path sorted_bam_index
        path genome
        path organelle
        path assembly_sizes
        val bam_bin_params
    output:
        path "output/*"
    script:
        output = "output"
    """
    #!/usr/bin/env bash

    mkdir -p $output
    echo "bam file : ${sorted_bam}, index file : ${sorted_bam_index}"
    s=\$(basename ${sorted_bam})
    regex="^(.+)-([^-]+)-Aligned[.]out[.]Sorted[.]bam\$"
    if [[ \$s =~ \$regex ]]; then
        assembly="\${BASH_REMATCH[1]}"
        run="\${BASH_REMATCH[2]}"
    else
        echo "Malformed BAM name, ${sorted_bam}"
        exit 1
    fi

    source $assembly_sizes
    total_size=\${assembly_sizes[\$assembly]}
    echo "Assembly \$assembly, run \$run, total size \$total_size"
    
    echo "genome ${genome}"
    head -5 ${genome}
    echo "organelle ${organelle}"
    head -5 ${organelle}
    echo "${genome}" > genome.mft
    echo "${organelle.join('\n')}" > organelle.mft
    samtools=`which samtools`
    bam_bin $bam_bin_params -bam $sorted_bam -o $output/\$assembly-\$run.bins -total-bam-size \$total_size \
        -fasta-manifest genome.mft -organelle-manifest organelle.mft -samtools-path \$samtools
    """

    stub:
        output = "output"
    """
    mkdir -p $output
    
    s=\$(basename ${sorted_bam})
    regex="^(.+)-([^-]+)-Aligned[.]out[.]Sorted[.]bam\$"
    if [[ \$s =~ \$regex ]]; then
        assembly="\${BASH_REMATCH[1]}"
        run="\${BASH_REMATCH[2]}"
    else
        echo "Malformed BAM name, ${sorted_bam}"
        exit 1
    fi

    touch $output/\$assembly-\$run.bins
    """
}


process merge_prepare {
    label 'large_disk'
    input:
        path runs
    output:
        path "*.bin_args"
    script:
    """
    #!/usr/bin/env bash

    declare -A bin_map
    for run in $runs; do
        s=\$(basename \${run})
        regex="^(.+)-([^-]+)[.]bins\$"
        if [[ \$s =~ \$regex ]]; then
            assembly="\${BASH_REMATCH[1]}"
        else
            echo "Malformed bins name, \${run}"
            exit 1
        fi
        for bam in \$run/*.bam; do
            key="\$assembly-\$(basename \$bam)"
            if [[ -z \${bin_map[\$key]} ]]; then
                bin_map[\$key]=\$bam
            else
                bin_map[\$key]="\${bin_map[\$key]} \$bam"
            fi
        done
    done
    for bin in \${!bin_map[@]}; do
        echo "merge --threads 8 \$bin \${bin_map[\$bin]}" >\$bin.bin_args
    done
    """

    stub:
    """
    touch 1.bin_args
    """
}


process merge {
    label 'big_job'
    label 'large_disk'
    input:
        path merge_args
        path bins
    output:
        path "*.bam"
    script:
    """
        samtools `cat $merge_args`   
    """
    
    stub:
    """
        touch 1.bam
    """
}
