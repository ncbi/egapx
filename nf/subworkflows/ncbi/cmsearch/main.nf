#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ----------------------------------------------------------------------------
workflow cmsearch_plane
{
take:
    genome_fasta  // path to genome fasta (can be gzipped)
                  // TODO: take asn_cache produced by upstream by setup_genome instead

    /*
    Expected params content:

    params:
        input:
            cmsearch: # these three files (can be local files or URLs)
              - https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/cmsearch/1/Rfam.seed
              - https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/cmsearch/1/rfam1410.cm
              - https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/support_data/cmsearch/1/rfam1410_amendments.xml
        tasks:
            cmsearch:
                num_batches: 16  # how many batches to split the genome into
                cmsearch_wnode: -exclusive-threshold 100000
                annot_merge: "-euk -unique"
    */

main:    
    prime_cache(genome_fasta)

    // Parallelized into num_batches - can be 1 if running locally;
    // otherwise should be approx chromosome-sized sequences, e.g. 16 or 32,
    // such that long-running jobs are split among different hosts.
    def num_batches = params.tasks?.cmsearch?.num_batches ?: 1

    cmsearch_results = cmsearch_wnode(
        prime_cache.out.asn_cache,
        prime_cache.out.seqids,
        params.input.cmsearch.files,
        batch_id=1..num_batches,
        num_batches=num_batches
    ) | collect
    
    gpx_qdump_and_annot_merge(prime_cache.out.asn_cache, cmsearch_results)

emit:
    cmsearch_annots = gpx_qdump_and_annot_merge.out.cmsearch_annots
}


// ----------------------------------------------------------------------------
// TODO: Consume asn_cache produced by upstream by setup_genome
process prime_cache
{
input:
    path "inp/fasta.fa"

output:
    path "out/asn_cache",    emit: asn_cache
    path "out/seqids.tsv",   emit: seqids

script:
    """
    set -exuo pipefail

    mkdir -p out/
    zcat -f inp/fasta.fa | prime_cache -ifmt fasta -cache out/asn_cache -oseq-ids out/seqids.tsv
    """
}

// ----------------------------------------------------------------------------
process cmsearch_wnode
{
    tag "${batch_id}" // NB: can only contain a number, otherwise egapx.py will crash in collect_logs()
                      // by incorrectly parsing task-name string that ends with a tag suffix.
    
input:
    path "inp/asn_cache"
    path "inp/seqids.tsv"
    
    // Passing these as inputs instead of accessing them directly from global params
    // because the paths may be URLs, and so they need to be staged by nextflow.
    path "inp/cmsearch_data/*"

    // invoknig cmsearch_node on each element of passed list [1..num_batches]
    each batch_id
    val num_batches
    
output:
    path "out/${batch_id}.gpx-job.asnb"

script:
    // TODO: is `which cmsearch` correct way to specify the path in egapx context,
    // or do we need an explicit path?
    """
    mkdir -p out

    #temporarily for development for faster turn-around
    #cat inp/seqids.tsv | tail -n 2 > tmp_seqids.tsv

    run_wnode_batch.py                                                  \\
        --exclusive                                                     \\
        --batch-num=${batch_id}                                         \\
        --num-batches=${num_batches}                                    \\
        --ids=inp/seqids.tsv                                            \\
        --asn-cache=inp/asn_cache                                       \\
        --work-dir=./var                                                \\
        --out-file=out/${batch_id}.gpx-job.asnb                         \\
        cmsearch_wnode                                                  \\
            -cmsearch-path   \$(dirname \$(which cmsearch))             \\
            -model-path      inp/cmsearch_data/rfam1410.cm              \\
            -rfam-amendments inp/cmsearch_data/rfam1410_amendments.xml  \\
            -rfam-stockholm  inp/cmsearch_data/Rfam.seed                \\
            -rfam-version    14.10                                      \\
            ${params.tasks?.cmsearch?.cmsearch_wnode ?: ''}             \\
    """
/*
    Notes on parallelization aspects of cmsearch:
    
  - cmsearch is the underlyig binary invoked by cmsearch_wnode.
    It is multithreaded (-cmsearch-cpu), but multithreading has 
    startup-overhead, and scales up to about 32 CPUs.

  - cmsearch_wnode invokes multiple cmsearch processes, controlled by
    -workers -workers-per-cpu -cpus-per-worker.

  - Jobs vary by orders of magnitute in size, e.g. one batch may have
    one job that is hundreds of MBs in size, and tens of thousands of 
    jobs 1kb in size.

  - Problem: if we use 1 CPU per cmsearch with many workers,
    then the runtime is bounded by the huge jobs tha cmseach could process
    faster if it ran in multithreaded mode. On the other extreme, a single
    multithreaded worker using all CPUs is bounded by the large number of
    small jobs, incurring the job-startup-overhead sequentially.

    The worker_node has ben modified to run in single-threaded-mode on short
    jobs and in multithreaded mode with 32 CPUs on large jobs exceeding
    100kb (-exclusive-threshold 100000). So we omit the -workers and 
    -cmsearch-cpu params and let the cmsearch_wnode manage it.

  - Problem: When running in cluster environment we want to split the jobs
    into equal-ish sized batches to be processed in paralell on multiple
    hosts, but we don't want to be processing them in parallel on a single
    host when running locally as that will overload the host.

    This can be managed in multiple ways:
        - the caller should specify num_batches=1 when running locally
        - the caller should specify maxForks=1 for cmsearch_wnode process.
          Note that this setting applies globally, whereas we are interested
          in maxForks=1 per-host (i.e. at most one cmsearch_wnode process
          per host).
        - run_wnode_batch.py --exclusive will invoke cmsearch_wnode 
          with flock, which will enforce sigle-process-per-host regardless.
*/

stub:
    """
    which run_wnode_batch.py
    mkdir -p ./out
    touch ./out/${batch_id}.gpx-job.asnb
    """
}


// ----------------------------------------------------------------------------
process gpx_qdump_and_annot_merge
{
input:
    path "inp/asn_cache"
    path 'inp/gpx/*'

output:
    path "out/cmsearch.asnb", emit: cmsearch_annots

script:
    """
    set -exuo pipefail

    mkdir -p out
    gpx_qdump -input-path ./inp/gpx/ -sort-by job-id -unzip '*' |
        annot_merge -asn-cache inp/asn_cache -output out/cmsearch.asnb ${params.tasks?.cmsearch?.annot_merge ?: ''}
    """

stub:
    """
    mkdir -p out
    touch out/cmsearch.asnb
    """
}
