#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow make_winmask_stats {
    take:
        genome_asnb
        seqids 
        gencoll_asn
        parameters
    main:
        stats = run_winmask_stats(genome_asnb, seqids, gencoll_asn, "").collect()
    emit:
        winmask_stats = stats
}


process run_winmask_stats {
    input:
        path genome_asnb
        path seqids
        path gencoll_asn
        val params
    output:
        path "output/*"
    script:
    genome_basename = genome_asnb.getBaseName().toString().replaceFirst(/\.(fa(sta)?|fna|asn(t|b)?)(\.gz)?$/, "")
    """
    echo $seqids > seqids.mft
    mkdir -p asncache
    prime_cache -cache asncache -ifmt asnb-seq-entry -i ${genome_asnb} -oseq-ids spids -split-sequences
    mkdir -p output
    make_winmask_stats -logfile ./mws.log -nogenbank -asn-cache asncache -gc-assembly ${gencoll_asn} -infmt seqids -input-manifest seqids.mft -out output/${genome_basename}.winmask_stats -sformat obinary -smem 1024
    cat ./mws.log
    """
    stub:
    genome_basename = genome_asnb.getBaseName().toString().replaceFirst(/\.(fa(sta)?|fna|asn(t|b)?)(\.gz)?$/, "")
    """
    mkdir -p output
    touch output/${genome_basename}.winmask_stats
    """
}
