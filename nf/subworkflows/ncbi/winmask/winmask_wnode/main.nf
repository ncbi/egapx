#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow winmask_wnode {
    take:
        genome_asnb
        seqids 
        gencoll_asn
        winmask_stats  
        parameters
    main:
        String gpx_qsubmit_params = merge_params("", parameters, "gpx_qsubmit")
        String winmasker_wnode_params = merge_params("", parameters, "winmasker_wnode")
        String gpx_make_outputs_params = merge_params("-unzip '*'", parameters, "gpx_make_outputs")
        String combine_blast_db_params = merge_params("", parameters, "combine_blast_db")
        
        (jobs, lines_per_file) = gpx_qsubmit(genome_asnb, seqids, gpx_qsubmit_params)
        jobs = jobs.flatten()
        wnode_collected = run_winmask_wnode(genome_asnb.collect(), winmask_stats.collect(), jobs, lines_per_file, winmasker_wnode_params).collect()
        mask = gpx_make_outputs(wnode_collected, gpx_make_outputs_params).collect()
        combine_blast_db(mask, combine_blast_db_params)
    emit:
        mask = combine_blast_db.out.mask
}


process gpx_qsubmit {
    label 'gpx_submitter'
    label 'small_mem'
    input:
        path genome_asnb
        path seqids
        val parameters
    output:
        path "job.*"
        env lines_per_file
    script:
        njobs=task.ext.split_jobs
    """
    echo $seqids > seqids.mft
    mkdir -p tmp/asncache
    prime_cache -cache tmp/asncache/ -ifmt asnb-seq-entry  -i ${genome_asnb} -oseq-ids spids -split-sequences

    gpx_qsubmit $parameters -ids-manifest seqids.mft -o jobs -nogenbank -asn-cache tmp/asncache/
    total_lines=\$(wc -l <jobs)
    (( lines_per_file = (\$total_lines + $njobs - 1) / $njobs ))
    echo total_lines=\$total_lines, lines_per_file=\$lines_per_file
    # split -l\$lines_per_file jobs job. -da 3
    # Use round robin to distribute jobs across nodes more evenly
    if [ \$total_lines -lt $njobs ]; then
        effective_njobs=\$total_lines
    else
        effective_njobs=$njobs
    fi
    split -nr/\$effective_njobs jobs job. -da 3
    rm -rf tmp
    """
    stub:
        njobs=16
    """
    for i in {1..$njobs}; do
        echo j.\${i} >> jobs
    done
    split -nr/$njobs jobs job. -da 3
    lines_per_file=10
    """
}



process run_winmask_wnode {
    label 'multi_node'
    label 'small_mem'
    input:
        path genome_asnb
        path winmask_stats // 
        path jobs
        val lines_per_file
        val parameters
    output:
        path "mask/*", emit: "mask"
    script:
    """
    mkdir -p tmp/interim
    mkdir -p tmp/asncache
    prime_cache -cache tmp/asncache/ -ifmt asnb-seq-entry  -i ${genome_asnb} -oseq-ids spids -split-sequences
    filename=\$(basename -- "$jobs")
    extension="\${filename##*.}"
    (( start_job_id = ((10#\$extension) * $lines_per_file) + 1 ))
    winmasker_wnode -ustat $winmask_stats -asn-cache tmp/asncache/ -workers ${task.ext.threads} -start-job-id \$start_job_id -input-jobs $jobs -nogenbank  -O tmp/interim $parameters > /dev/null 2> /dev/null
    mkdir -p mask
    cat tmp/interim/* > mask/winmasker_wnode.${task.index}.gpx-job.asnb
    rm -rf tmp
    """
    stub:
    """
    mkdir -p mask
    touch mask/winmasker_wnode.${task.index}.gpx-job.asnb
    """
}


process gpx_make_outputs {
    label 'single_cpu'
    label 'small_mem'
    input:
        path files, stageAs: "gpx_inputs/*"
        val parameters
    output:
        path "output/*", emit: "all"
    script:
    """
    ls -1 gpx_inputs/* > gpx_inputs.mft
    mkdir -p output
    gpx_make_outputs $parameters  -input-manifest gpx_inputs.mft -output output/winmask.concat.#.asnb -output-manifest winmask.concat.mft  -num-partitions 1
    """
    stub:
    """
    mkdir -p output
    touch output/winmask.asnb
    """
}


process combine_blast_db {
    label 'single_cpu'
    label 'small_mem'
    input:
        path files, stageAs: "combine_db_inputs/*"
        val parameters
    output:
        path "winmask.asnb", emit: "mask"
    script:
    """
    ls -1 combine_db_inputs/* > winmask.concat.mft
    combine_blast_db -input-manifest winmask.concat.mft -o winmask.asnb
    """
    stub:
    """
    touch winmask.asnb
    """
}
