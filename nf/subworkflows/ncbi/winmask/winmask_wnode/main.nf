#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'

split_count=16


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
    input:
        path genome_asnb
        path seqids
        val params
    output:
        path "job.*"
        env lines_per_file
    script:
        njobs=split_count
    """
    echo $seqids > seqids.mft
    mkdir -p ./asncache/
    prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ${genome_asnb} -oseq-ids spids -split-sequences

    gpx_qsubmit $params -ids-manifest seqids.mft -o jobs -nogenbank -asn-cache ./asncache/ 
    total_lines=\$(wc -l <jobs)
    (( lines_per_file = (total_lines + ${njobs} - 1) / ${njobs} ))
    echo total_lines=\$total_lines, lines_per_file=\$lines_per_file
    # split -l\$lines_per_file jobs job. -da 3
    # Use round robin to distribute jobs across nodes more evenly
    if [ \$total_lines -lt $njobs ]; then
        effective_njobs=\$total_lines
    else
        effective_njobs=$njobs
    fi
    split -nr/\$effective_njobs jobs job. -da 3
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
    input:
        path genome_asnb
        path winmask_stats // 
        path jobs
        val lines_per_file
        val parameters
    output:
        path "mask/*", emit: "mask"
    script:
        job_num = jobs.toString().tokenize('.').last().toInteger()
    """
    njobs=`wc -l <$jobs`
    if [ \$njobs -lt 16 ]; then
        threads=\$njobs
    else
        threads=16
    fi
    mkdir -p interim
    mkdir -p asncache
    prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ${genome_asnb} -oseq-ids spids -split-sequences
    filename=\$(basename -- "$jobs")
    extension="\${filename##*.}"
    (( start_job_id = ((10#\$extension) * $lines_per_file) + 1 ))
    winmasker_wnode -ustat $winmask_stats -asn-cache ./asncache/ -workers \$threads -start-job-id \$start_job_id -input-jobs $jobs -nogenbank  -O interim $parameters
    mkdir -p mask
    for f in interim/*; do
        if [ -f \$f ]; then
            mv \$f mask/\${extension}_\$(basename \$f)
        fi
    done
    """
    stub:
        job_num = jobs.toString().tokenize('.').last().toInteger()
    """
    mkdir -p mask
    touch mask/${job_num}.asnb
    """
}


process gpx_make_outputs {
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
