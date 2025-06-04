#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'

// Get protein FASTA from final annotation

workflow generate_fasta_from_annots {
    take:
        gencoll_asn
        scaffolds    // asn seqentrs
        chromosomes  // asn seqentrysseqids
        parameters   // Map : extra parameter and parameter update
    main:
        gnomon_asn2fasta(gencoll_asn, scaffolds, chromosomes, merge_params('', parameters, 'gnomon_asn2fasta'))
    emit:
        proteins = gnomon_asn2fasta.out.proteins
}


process gnomon_asn2fasta {
    input:
        path gencoll_asn, stageAs: 'gencoll.asn'
        path scaffolds, stageAs: 'input/*' // asn seqentry
        path chromosomes, stageAs: 'input/*' // asn seqentry
        val params
    output:
        path 'output/*.prot.fa', emit: 'proteins'
    script:
    """
    mkdir -p output
    mkdir -p asncache
    for f in input/*; do
        echo \$f >> input.mft
    done
    prime_cache -cache asncache -input-manifest input.mft -ifmt asn-seq-entry
    gnomon_asn2fasta -nogenbank -asn-cache asncache -gc-assembly gencoll.asn -ifmt seq-entry -input-manifest input.mft -prot-output output/@.prot.fa -prot-output-manifest output/proteins.mft $params
    """
    stub:
    """
    mkdir -p output
    touch output/assembly.prot.fa
    """
}