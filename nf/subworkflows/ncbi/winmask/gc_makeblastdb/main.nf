#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'

shrd_mft_name="blastdb.mft"    // this must be shared with tblastn_align

workflow gc_makeblastdb {
    take:
        genome_asnb
        seqids 
        gencoll_asn
        alternate_softmask
        parameters
    main:
        blastdb = run_gc_makeblastdb(genome_asnb, seqids, gencoll_asn, alternate_softmask, "")
    emit:
        blastdb = blastdb
}


process run_gc_makeblastdb{
    label 'single_cpu'
    label 'small_mem'
    input:
        path genome_asnb
        path seqids
        path gencoll_asn
        path alternate_softmask
        val params
    output:
        path "output/*", emit: 'blastdb'   // this is not the output directory itself, it is the collection of files under output. 
                                           // process that use this should expect to get the collection of blastdb files staged as they direct
                                           // without another directory layer around it. 
    script:
    """
    echo "${alternate_softmask.join('\n')}" > softmask_data.mft
    mkdir -p tmp/asncache
    prime_cache -cache tmp/asncache -ifmt asnb-seq-entry -i ${genome_asnb} -oseq-ids spids -split-sequences
    mkdir -p output
    gc_makeblastdb -nogenbank -asn-cache tmp/asncache -gc-assembly $gencoll_asn -input $seqids -softmask-manifest softmask_data.mft -output-manifest output/$shrd_mft_name -output-path output -title 'BLASTdb created by EGAPx'
    rm -rf tmp
    """
    stub:
    """
    mkdir -p output
    touch output/blastdb.nal
    echo "blastdb" > output/blastdb.mft 
    """
}
