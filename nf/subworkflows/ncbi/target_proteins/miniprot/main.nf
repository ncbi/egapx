#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow miniprot {
    take:
        fasta_genome_file  //path: genome fasta file
        fasta_proteins_file  //path: protein fasta file
        parameters      // Map : extra parameter and parameter update
    main:
        default_params = "-t 8"
        def value = parameters.get('miniprot', "")
        value = value.replaceFirst("-cpu-count", "-t")
        value = value.replaceFirst("-max-intron", "-G")
        parameters['miniprot'] = value
        effective_params = merge_params(default_params, parameters, 'miniprot')
        run_miniprot(fasta_genome_file, fasta_proteins_file, effective_params)

    emit:
        miniprot_file = run_miniprot.out.miniprot_file
}


process run_miniprot {
    input:
        path fasta_genome_file
        path fasta_proteins_file
        val parameters
    output:
        path ('output/aligns.paf'), emit: 'miniprot_file'

    script:
    """
    mkdir -p output
    miniprot ${parameters}  ${fasta_genome_file} ${fasta_proteins_file} > output/aligns.paf
    """
}
