#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'

bdb_dir = "blast__db"
shrd_mft_name="blastdb.mft"    // this must be shared with gc_makeblastdb
mft_name = "$bdb_dir/$shrd_mft_name"

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
    label 'gpx_submitter'
    label 'small_mem'
    input:
        path genome_asn
        path proteins_asn
        path blastdb, stageAs: bdb_dir+"/*"
        val params
    output:
        path "job.*"
        env lines_per_file
    script:
        njobs=task.ext.split_jobs
    """
    mkdir -p tmp/asncache
    auto_prime_cache.py -cache tmp/asncache/ -i ${proteins_asn} -oseq-ids ./prot_ids -split-sequences

    # NB: Use db alias file, *.nal without extension to set the database. 
    # Use -db instead of -db-manifest, the latter generates absolute paths 
    # for the database and will not work on clusters without shared filesystem    
    #gpx_qsubmit $params -ids prot_ids -o jobs -nogenbank -asn-cache tmp/asncache/ -db blastdb/blastdb/blastdb
    #gpx_qsubmit $params -ids prot_ids -o jobs -nogenbank -asn-cache tmp/asncache/ -db-manifest blastdb/blastdb.mft
    db_name=\$( cat $mft_name | grep  -o '[^/]*\$' )
    db_path_name="$bdb_dir/\${db_name}"
    gpx_qsubmit $params -ids prot_ids -o jobs -nogenbank -asn-cache tmp/asncache/ -db \${db_path_name}
     
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


process run_tblastn_wnode {  //TODO wnode 
    label 'long_job'
    label 'multi_node'
    label 'med_mem'
    input:
        path genome_asn
        path proteins_asn
        path blastdb, stageAs: bdb_dir+"/*"
        path jobs
        val lines_per_file
        val parameters
    output:
        path "blast/*", emit: "blast"
    script:
    """
    mkdir -p tmp/interim
    mkdir -p tmp/asncache

    auto_prime_cache.py -cache tmp/asncache/ -i ${genome_asn} -oseq-ids spids -split-sequences
    auto_prime_cache.py -cache tmp/asncache/ -i ${proteins_asn} -oseq-ids spids2 -split-sequences
    filename=\$(basename -- "$jobs")
    extension="\${filename##*.}"
    (( start_job_id = ((10#\$extension) * $lines_per_file) + 1 ))
    tblastn_wnode -asn-cache tmp/asncache/ -workers ${task.ext.threads} -start-job-id \$start_job_id -input-jobs $jobs -nogenbank -O tmp/interim $parameters
    
    mkdir -p blast
    cat tmp/interim/* > blast/tblastn_wnode.${task.index}.gpx-job.asnb
    rm -rf tmp
    """
    stub:
    """
    mkdir -p blast
    touch blast/tblastn_wnode.${task.index}.gpx-job.asnb
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
    gpx_make_outputs $parameters  -input-manifest gpx_inputs.mft -output output/blast.#.asn -output-manifest blast_align.mft  -num-partitions 1
    """
    stub:
    """
    mkdir -p output
    touch output/blast.asn
    """
}
