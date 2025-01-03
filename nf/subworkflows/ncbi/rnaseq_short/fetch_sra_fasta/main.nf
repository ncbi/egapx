#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'
params.inpdir = ""
workflow fetch_sra_fasta  {
    take:
        sra_run_file  // Channel: sra run file
        parameters     // Map : extra parameter and parameter update
    main:
        a = read_sra_file(sra_run_file)
        b =  a.map { it -> it.split()}
        c =  run_fetch_sra_fasta(b.flatten())

    emit:
        //fasta_pair_list = run_fetch_sra_fasta.out.fasta_pair_list
        fasta_pair_list = c
}


process read_sra_file {
    input:
        path sra_run_file
    output:
        env  exitvar
    script:

    """
    exitvar=()
    while read -r line; do  [[ \$line = \\#* ]] && continue; exitvar+=(\"\$line\"); done < ${sra_run_file}
    """

    stub:
    """
    exitvar=SRA000001
    """
}


process run_fetch_sra_fasta {
    input:
        val sra
    output:
        tuple val (sra),  path ('output/*.{1,2}')  , emit: 'fasta_pair_list'
    script:
    """
    output_${sra}=\$(srapath ${sra})
    curl -o ${sra} \$output_${sra}
    fasterq-dump --skip-technical --threads 6 --split-files --seq-defline ">\\\$ac.\\\$si.\\\$ri" --fasta -O .  ./${sra}
    rm -f ${sra}
    mkdir -p output
    if [ -f ${sra}_1.fasta ]; then
        mv ${sra}_1.fasta output/${sra}.1
    fi
    if [ -f ${sra}_2.fasta ]; then
        mv ${sra}_2.fasta output/${sra}.2
    fi
    """
    stub:
    """
    mkdir -p output
    touch output/${sra}.1
    touch output/${sra}.2
    """
}
