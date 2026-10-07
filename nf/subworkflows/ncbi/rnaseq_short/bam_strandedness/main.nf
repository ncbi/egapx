#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow bam_strandedness {
    take:
        bam_list              // list: BAM
        sra_metadata        // path: file with sra metadata
        parameters          // Map : extra parameter and parameter update.
    main:
        // Not yet used but supposed to be used by rnaseq_divide_by_strandedness
        rnaseq_divide_by_strandedness_params = merge_params("-min-aligned 1000000 -min-unambiguous 200 -min-unambiguous-pct 2 -max-unambiguous-pct 100 -percentage-threshold 98", parameters, 'rnaseq_divide_by_strandedness')
        rnaseq_divide_by_strandedness(bam_list, sra_metadata, rnaseq_divide_by_strandedness_params)
    emit:
        strandedness = rnaseq_divide_by_strandedness.out.strandedness
        normalized_strandedness = rnaseq_divide_by_strandedness.out.normalized_strandedness
        stranded_runs = rnaseq_divide_by_strandedness.out.normalized_stranded_runs
        unstranded_runs = rnaseq_divide_by_strandedness.out.normalized_unstranded_runs
        all = rnaseq_divide_by_strandedness.out.all
}


process rnaseq_divide_by_strandedness {
    label 'large_disk'
    label 'single_cpu'
    label 'small_mem'
    input:
        path bam_list
        path metadata_file
        val parameters
    output:
        path "output/run.strandedness", emit: 'strandedness'
        path "normalized/run.strandedness", emit: 'normalized_strandedness'
        path "output/stranded.list", emit: 'stranded_runs', optional: true
        path "output/unstranded.list", emit: 'unstranded_runs', optional: true
        path "normalized/stranded.list", emit: 'normalized_stranded_runs', optional: true
        path "normalized/unstranded.list", emit: 'normalized_unstranded_runs', optional: true
        path "output/*", emit: 'all' 
    script:
    """
    mkdir -p output
    mkdir -p normalized
    mkdir -p tmp
    samtools=\$(which samtools)
    if [ $bam_list == "unpacked_genome.bam" ]; then
        mv unpacked_genome.bam GCF_030936135.1_lcl-SRR10853086-Aligned.out.bam
        echo "GCF_030936135.1_lcl-SRR10853086-Aligned.out.bam" > bam_list.mft
    else
        echo "${bam_list.join('\n')}" > bam_list.mft
    fi
    rnaseq_divide_by_strandedness -work-area tmp -align-manifest bam_list.mft -metadata $metadata_file  $parameters  -samtools-executable \$samtools -stranded-output output/stranded.list -strandedness-output output/run.strandedness -unstranded-output output/unstranded.list

    normalize_run_ids() {
        local input_file="\$1"
        local output_file="\$2"
        if [ -f "\$input_file" ]; then
            awk -F '\\t' 'BEGIN { OFS = "\\t" }
                FNR == NR {
                    if (\$0 !~ /^#/ && \$1 != "") known[\$1] = 1
                    next
                }
                \$0 ~ /^#/ { print; next }
                {
                    run_id = \$1
                    if (!(run_id in known)) {
                        best = ""
                        for (accession in known) {
                            if (index(run_id, accession "_") == 1 && length(accession) > length(best)) {
                                best = accession
                            }
                        }
                        if (best == "") {
                            printf "Cannot map strandedness run ID %s to an SRA accession in %s\\n", run_id, FILENAME > "/dev/stderr"
                            exit 1
                        }
                        \$1 = best
                    }
                    print
                }
            ' "$metadata_file" "\$input_file" > "\$output_file"
        fi
    }

    normalize_run_ids output/run.strandedness normalized/run.strandedness
    normalize_run_ids output/stranded.list normalized/stranded.list
    normalize_run_ids output/unstranded.list normalized/unstranded.list

    rm -rf tmp
    """
    stub:
    """
    mkdir -p output
    mkdir -p normalized
    touch output/run.strandedness
    touch output/stranded.list
    touch output/unstranded.list
    touch normalized/run.strandedness
    touch normalized/stranded.list
    touch normalized/unstranded.list
    """
}
