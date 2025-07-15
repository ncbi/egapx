#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'

split_count=4


workflow tblastn_align {
    take:
        genome_asn
        proteins_asn
        blastdb
        parameters
    main:
        String gpx_qsubmit_params = merge_params("", parameters, "gpx_qsubmit")
        String tblastn_wnodeparams = merge_params("", parameters, "tblastn_wnode")
        String gpx_make_outputs_params = merge_params("-unzip '*'", parameters, "gpx_make_outputs")
        
        (jobs, lines_per_file) = gpx_qsubmit(genome_asn, proteins_asn, blastdb, gpx_qsubmit_params)
        jobs = jobs.flatten()
        wnode_collected = run_tblastn_wnode(genome_asn.collect(), proteins_asn.collect(),  blastdb.collect(), jobs, lines_per_file, tblastn_wnodeparams).collect()
        blast_asn = gpx_make_outputs(wnode_collected, gpx_make_outputs_params).collect()
    emit:
        blast_asn = blast_asn
}


process gpx_qsubmit {
    input:
        path genome_asn
        path proteins_asn
        path blastdb
        val params
    output:
        path "job.*"
        env lines_per_file
    script:
        njobs=split_count
    """
    echo 'prot_ids' > gilist.mft
    echo $blastdb/blastdb > blastdb.mft
    
    mkdir -p ./asncache/
    prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ${proteins_asn} -oseq-ids ./prot_ids -split-sequences
    
    gpx_qsubmit $params -ids-manifest gilist.mft -o jobs -nogenbank -asn-cache ./asncache/  -db-manifest blastdb.mft   
    
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



process run_tblastn_wnode {
    input:
        path genome_asn
        path proteins_asn
        path blastdb  
        path jobs
        val lines_per_file
        val parameters
    output:
        path "blast/*", emit: "blast"
    script:
    """
    njobs=`wc -l <$jobs`
    if [ \$njobs -lt 16 ]; then
        threads=\$njobs
    else
        threads=4
    fi
    threads=8
    mkdir -p interim
    mkdir -p ./asncache/
    prime_cache -cache ./asncache/ -ifmt asn-seq-entry  -i ${genome_asn} -oseq-ids spids -split-sequences
    prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ${proteins_asn} -oseq-ids spids2 -split-sequences
    filename=\$(basename -- "$jobs")
    extension="\${filename##*.}"
    (( start_job_id = ((10#\$extension) * $lines_per_file) + 1 ))
    tblastn_wnode -asn-cache  ./asncache/   -workers \$threads -start-job-id \$start_job_id -input-jobs $jobs -nogenbank  -O interim $parameters 
    
    mkdir -p blast
    cat interim/* > blast/tblastn_wnode.${task.index}.gpx-job.asnb
    rm -rf interim
    """
    stub:
    """
    mkdir -p blast
    touch blast/tblastn_wnode.${task.index}.gpx-job.asnb
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
    gpx_make_outputs $parameters  -input-manifest gpx_inputs.mft -output output/blast.#.asn -output-manifest blast_align.mft  -num-partitions 1
    """
    stub:
    """
    mkdir -p output
    touch output/blast.asn
    """
}

