#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../utilities'

params.verbose = false

// To save staging of the lineage files we pass only the specific lineage
// and arrange staging so that it results in an order expected by BUSCO, i.e.
// the busco_download directory contains this specific lineage in the subdirectory
// with the lineage name
workflow busco {
    take: 
        proteins         // path: FASTA with proteins
        lineage          // val: lineage name
        lineage_download // path: download directory with the lineage
        parameters       // Map : extra parameter and parameter update
    main:
        run_busco(proteins, lineage, lineage_download, parameters)
    emit:
        results = run_busco.out.results
}


process run_busco {
    input:
        path proteins
        val  lineage
        path lineage_download
        val  parameters
    output:
        path "output/*", emit: 'results'
    script:
        // For some executors (e.g. SGE) the task.cpus is not set correctly, they allocate correct number of threads through clusterOptions.
        // We use ext.cpus to pass the number of cpus here and use clusterOptions to allocate large enough instance
        def cpus = task.cpus == 1 && task.ext.cpus ? task.ext.cpus : task.cpus
    """
    download_params=''
    if [ -d "$lineage_download" ]; then
        mkdir -p busco_downloads/lineages
        ln -s `readlink -f $lineage_download` busco_downloads/lineages/$lineage
        download_params="--download_path busco_downloads --offline"
    fi
    if ! busco \$download_params -i $proteins --out output --mode proteins --cpu $cpus --lineage_dataset $lineage; then
        mkdir -p output
        cp busco*.log output | true
        touch output/short_summary.json
        touch output/short_summary.txt
    fi
    """
    stub:
    """
    mkdir -p output
    touch output/short_summary.json
    touch output/short_summary.txt
    """
}
