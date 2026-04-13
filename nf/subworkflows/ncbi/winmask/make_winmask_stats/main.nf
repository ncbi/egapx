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
        String make_stats_params = merge_params(" -infmt seqids -sformat binary -mem 1024 -smem 1024 -use-stored-stats never", parameters, "make_winmask_stats")
        stats = run_winmask_stats(genome_asnb, seqids, gencoll_asn, make_stats_params).collect()
    emit:
        winmask_stats = stats
}


process run_winmask_stats {
    label 'single_cpu'
    label 'small_mem'
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
    mkdir -p tmp/asncache
    prime_cache -cache tmp/asncache -ifmt asnb-seq-entry -i ${genome_asnb} -oseq-ids spids -split-sequences
    mkdir -p output
    make_winmask_stats -logfile ./mws.log -nogenbank -asn-cache tmp/asncache -gc-assembly ${gencoll_asn} -input-manifest seqids.mft -out output/${genome_basename}.winmask_stats  ${params}
    cat ./mws.log
    rm -rf tmp
    """
    stub:
    genome_basename = genome_asnb.getBaseName().toString().replaceFirst(/\.(fa(sta)?|fna|asn(t|b)?)(\.gz)?$/, "")
    """
    mkdir -p output
    touch output/${genome_basename}.winmask_stats
    """
}
