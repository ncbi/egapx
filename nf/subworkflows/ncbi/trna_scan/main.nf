#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ----------------------------------------------------------------------------
workflow trnascan_plane
{
take:
    genome_fasta  // path to genome fasta (can be gzipped)
                  // TODO: take asn_cache produced by upstream by setup_genome instead
main:    
    prime_cache(genome_fasta)

    def num_batches = params?.tasks?.trnascan?.num_batches ?: 1

    trnascan_results = trnascan_wnode(
        prime_cache.out.asn_cache,
        prime_cache.out.seqids,
        batch_id=1..num_batches,
        num_batches=num_batches
    ) | collect
    
    trnascan_dump(prime_cache.out.asn_cache, trnascan_results)

emit:
    trnascan_annots = trnascan_dump.out.trnascan_annots
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
process trnascan_wnode
{
    tag "${batch_id}" // NB: can only contain a number, otherwise egapx.py will crash in collect_logs()
                      // by incorrectly parsing task-name string that ends with a tag suffix.
    
    ext.trnascan_params = "-X 55 -Q -b -q"
input:
    path "inp/asn_cache"
    path "inp/seqids.tsv"
    
    // invoknig trnascan_node on each element of passed list [1..num_batches]
    each batch_id
    val num_batches
    
output:
    path "out/${batch_id}.gpx-job.asnb"

script:
    """
    mkdir -p out

    run_wnode_batch.py                                                  \\
        --exclusive                                                     \\
        --batch-num=${batch_id}                                         \\
        --num-batches=${num_batches}                                    \\
        --ids=inp/seqids.tsv                                            \\
        --asn-cache=inp/asn_cache                                       \\
        --work-dir=./var                                                \\
        --out-file=out/${batch_id}.gpx-job.asnb                         \\
        trnascan_wnode                                                  \\
            ${task.ext.trnascan_params} -tRNAscan \$(which tRNAscan-SE) \\
    """
stub:
    """
    mkdir -p ./out
    touch ./out/${batch_id}.gpx-job.asnb
    """
}


// ----------------------------------------------------------------------------
process trnascan_dump
{
    ext.trnascan_dump_params = "-X 55"

input:
    path "inp/asn_cache"
    path 'inp/gpx/*'

output:
    path "out/trnascan.asnb", emit: trnascan_annots

script:
    """
    set -exuo pipefail

    mkdir -p out var
    gpx_qdump -input-path ./inp/gpx/ -sort-by job-id -unzip '*' |
        trnascan_dump -oasn out/trnascan.asnb -ostruc var/struc.tar.gz ${task.ext.trnascan_dump_params}
    """

stub:
    """
    mkdir -p out
    touch out/trnascan.asnb
    """
}
