#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow prosplign_wnode {
    take:
        genome_asnb
        proteins_asnb 
        compartments
        max_intron  
        parameters
    main:
        String gpx_qsubmit_params = merge_params("", parameters, "gpx_qsubmit")
        String prosplign_wnode_params = merge_params("", parameters, "prosplign_wnode")
        String gpx_qdump_params = merge_params("-unzip '*'", parameters, "gpx_qdump")
        
        (jobs, lines_per_file) = gpx_qsubmit(genome_asnb, proteins_asnb, compartments, gpx_qsubmit_params)
        jobs = jobs.flatten()
        wnode_collected = run_prosplign_wnode(genome_asnb.collect(),  proteins_asnb.collect(), compartments.collect(), jobs, max_intron, lines_per_file.first(), prosplign_wnode_params).collect()
        prosplign = gpx_qdump(wnode_collected, gpx_qdump_params).collect()
    emit:
        prosplign_asn = prosplign
}


process gpx_qsubmit {
    label 'gpx_submitter'
    label 'small_mem'
    input:
        path genome_asnb
        path proteins_asnb
        path compartments
        val params
    output:
        path "job.*"
        env lines_per_file
    script:
        njobs=task.ext.split_jobs
    """
    mkdir -p tmp/asncache
    auto_prime_cache.py -cache tmp/asncache/  -i ${genome_asnb} -oseq-ids spids -split-sequences
    auto_prime_cache.py -cache tmp/asncache/ -i ${proteins_asnb} -oseq-ids spids2 -split-sequences
    gpx_qsubmit $params -asn $compartments -o jobs -nogenbank -asn-cache tmp/asncache/ 
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



process run_prosplign_wnode {
    label 'long_job'
    label 'multi_node'  //this is not determined yet. subject to change
    label 'med_mem'
    input:
        path genome_asnb
        path proteins_asnb
        path compartments
        path jobs
        val max_intron
        val lines_per_file
        val parameters
    output:
        path "output/*", emit: "prosplign"
    script:
        job_num = jobs.toString().tokenize('.').last().toInteger()
    """
    mkdir -p tmp/interim
    mkdir -p tmp/asncache
    auto_prime_cache.py -cache tmp/asncache/ -i ${genome_asnb} -oseq-ids spids -split-sequences
    auto_prime_cache.py -cache tmp/asncache/ -i ${proteins_asnb} -oseq-ids spids2 -split-sequences
    
    filename=\$(basename -- "$jobs")
    extension="\${filename##*.}"
    (( start_job_id = ((10#\$extension) * $lines_per_file) + 1 ))
    prosplign_wnode -max_intron $max_intron -asn-cache tmp/asncache/ -workers ${task.ext.threads} -start-job-id \$start_job_id -input-jobs $jobs -nogenbank  -O tmp/interim $parameters
    mkdir -p output
    cat tmp/interim/* > output/prosplign_wnode.${task.index}.gpx-job.asnb
    rm -rf tmp
    """

    stub:
    """
    mkdir -p output
    echo ${task.index} > output/prosplign_wnode.${task.index}.gpx-job.asnb
    """
}


process gpx_qdump {
    label 'single_cpu'
    label 'small_mem'
    input:
        path files, stageAs: "inputs/*"
        val params
    output:
        path "*.asn", emit: "prosplign_asn"
    script:
    """
    gpx_qdump $params  -input-path inputs -output prosplign.asn
    """
    stub:
    """
    touch prosplign.asn
    """
}
