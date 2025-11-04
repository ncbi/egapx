#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'

split_count=16


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
    input:
        path genome_asnb
        path proteins_asnb
        path compartments
        val params
    output:
        path "job.*"
        env lines_per_file
    script:
        njobs=split_count
    """
    mkdir -p ./asncache/
    prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ${genome_asnb} -oseq-ids spids -split-sequences
    prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ${proteins_asnb} -oseq-ids spids2 -split-sequences
    gpx_qsubmit $params -asn $compartments -o jobs -nogenbank -asn-cache ./asncache/ 
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



process run_prosplign_wnode {
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
    njobs=`wc -l <$jobs`
    if [ \$njobs -lt 16 ]; then
        threads=\$njobs
    else
        threads=16
    fi
    mkdir -p interim
    mkdir -p asncache
    prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ${genome_asnb} -oseq-ids spids -split-sequences
    prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ${proteins_asnb} -oseq-ids spids2 -split-sequences
    
    filename=\$(basename -- "$jobs")
    extension="\${filename##*.}"
    (( start_job_id = ((10#\$extension) * $lines_per_file) + 1 ))
    prosplign_wnode -max_intron $max_intron -asn-cache ./asncache/ -workers \$threads -start-job-id \$start_job_id -input-jobs $jobs -nogenbank  -O interim $parameters
    mkdir -p output
    cat interim/* > output/prosplign_wnode.${task.index}.gpx-job.asnb
    rm -rf interim
    """

    stub:
    """
    mkdir -p output
    echo ${task.index} > output/prosplign_wnode.${task.index}.gpx-job.asnb
    """
}


process gpx_qdump {
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






/*gpx_qsubmit -dryrun -NxM-threshold 10000 -affinity subject -asn-cache sequence_cache -asn-manifest compartments.mft
prosplign_wnode -asn-cache sequence_cache -backlog 1 -compartment_penalty 0.5 -cut_flank_partial_codons true -fill_holes false -flank_positives 55 
-frameshift_opening 30 -gap_extension 1 -gap_opening 10 -intron_AT 25 -intron_GC 20 -intron_GT 15 -intron_non_consensus 34 -inverted_intron_extension 1000
 -max_bad_len 45 -max_extent 10000 -max_intron 400000 -maximize coverage -min_compartment_idty 0.5 -min_exon_ident 30 -min_exon_positives 55 
 -min_flanking_exon_len 15 -min_good_len 59 -min_intron_len 30 -min_positives 15 -min_singleton_idty 0.65 -service GPipeExec_Prod -start_bonus 8 
 -stop_bonus 8 -total_positives 70 -output prosplign.asn -output-manifest prosplign_align.mft -unzip '*'  -destination-dir out/
gpx_qdump -output prosplign.asn -output-manifest prosplign_align.mft -unzip '*'  -destination-dir out/  */
