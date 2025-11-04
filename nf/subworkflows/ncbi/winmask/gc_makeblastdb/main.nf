#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


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
    input:
        path genome_asnb
        path seqids
        path gencoll_asn
        path alternate_softmask
        val params
    output:
        path "output/*", emit: 'blastdb'
    script:
    """
    echo $seqids > seqids.mft
    echo $gencoll_asn > gencoll_asn.mft
    echo "${alternate_softmask.join('\n')}" > softmask_data.mft
    mkdir -p asncache
    prime_cache -cache asncache -ifmt asnb-seq-entry -i ${genome_asnb} -oseq-ids spids -split-sequences
    mkdir -p output
    gc_makeblastdb -nogenbank -asn-cache asncache -gc-assembly-manifest gencoll_asn.mft -input-manifest seqids.mft -softmask-manifest softmask_data.mft -output-path output  -title 'BLASTdb created by EGapx' -output-assembly-unit-manifest output/blastdb_assm_unit.mft  -output-full-assembly-manifest output/blastdb_assm.mft  -output-manifest output/blastdb.mft 
    """
    stub:
    """
    mkdir -p output
    touch output/blastdb.nal
    """
}
