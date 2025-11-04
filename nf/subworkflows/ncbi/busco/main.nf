#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../utilities'

params.verbose = false

// To save staging of the lineage files we pass only the specific lineage
// and arrange staging so that it results in an order expected by BUSCO, i.e.
// the busco_download directory contains this specific lineage in the subdirectory
// with the lineage name.
workflow busco {
    take: 
        proteins         // path: FASTA with proteins.
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
    """
    download_params=''
    if [ -d "$lineage_download" ]; then
        mkdir -p busco_downloads/lineages
        ln -s `readlink -f $lineage_download` busco_downloads/lineages/$lineage
        download_params="--download_path busco_downloads --offline"
    fi

    num_threads=\$((\$(nproc) < 64 ? \$(nproc) : 64))

    # Set OPENBLAS_NUM_THREADS=1 because otherwise the number of threads is 32*num_threads,
    # which will hit the ulimit (typically 1024). That's on top of nextflow itself using ~200 threads.
    export OPENBLAS_NUM_THREADS=1

    if ! busco \$download_params -i $proteins --out output --mode proteins --cpu \$num_threads --lineage_dataset $lineage --tar; then

        # if busco errored-out, that's not a fatal failure for the pipeline.
        # Instead, just add stub outputs. NB: log, if exists, is in output/logs/busco.log
        mkdir -p output
        touch output/short_summary.json
        touch output/short_summary.txt
    fi

    # this contains tens of thousands of files, so clean it up.
    rm -rf busco_downloads 
    """

    stub:
    """
    mkdir -p output
    touch output/short_summary.json
    touch output/short_summary.txt
    """
}
