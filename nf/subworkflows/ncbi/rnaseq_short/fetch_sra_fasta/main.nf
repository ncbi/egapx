#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// TODO: Move to utilities or other common location

include { merge_params } from '../../utilities'
params.inpdir = ""
workflow fetch_sra_fasta  {
    take:
        sra_id_channel // Channel of SRA IDs (strings)
        parameters     // Map : parameter update

    main:
        def split_reads = parameters.get('split', '')

        // Pass the channel of IDs directly to run_fetch_sra_fasta
        c = run_fetch_sra_fasta(sra_id_channel, split_reads)

        paired = c.fasta_pair_list.flatMap { sra, _1_files, _2_files ->
            // ensure files are sorted so that corresponding parts line up
            def list_1 = _1_files instanceof List ? _1_files.sort() : ( _1_files ? [ _1_files ] : [])
            def list_2 = _2_files instanceof List ? _2_files.sort() : ( _2_files ? [ _2_files ] : [])
            if (list_2.isEmpty()) {
                // Single-end: just emit _1 files individually
                list_1.collect { f1 -> tuple(sra, [f1]) }
            } else {
                // Paired-end: pair up corresponding files
                [list_1, list_2].transpose().collect { pair -> 
                    tuple(sra, pair)
                }
            }
        }

    emit:
        fasta_pair_list = paired
}

process read_sra_file {
    label 'single_cpu'
    label 'small_mem'
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
    exitvar="SRA000001 SRA000002"
    """
}

process run_fetch_sra_fasta {
    label 'single_cpu'
    label 'small_mem'
    input:
        val sra
        val split_reads
    output:
        tuple val (sra), path ('output/*_1*.fasta*'), path ('output/*_2*.fasta*', arity: '0..*'), emit: 'fasta_pair_list'
    script:
    def split_cmd = split_reads ? """
    for f in output/${sra}_*.fasta; do
        seqkit split2 -p ${split_reads} -O output_split "\$f"
    done
    rm -f output/${sra}_*.fasta
    mv output_split/* output/
    rmdir output_split
    gzip output/*.fasta
    """ : ""
    """
    mkdir -p output
    curl -fL --retry 5 -C - -o ${sra}.sra \$(srapath ${sra})
    fasterq-dump --skip-technical --threads 6 --split-files --seq-defline ">gnl|SRA|\\\$ac.\\\$si.\\\$ri" --fasta --outdir output  ./${sra}.sra
    ls output/${sra}_*.fasta
    rm -f ${sra}.sra
    ${split_cmd}
    """
    stub:
    def stub_split_cmd = split_reads ? """
    for i in \$(seq 1 ${split_reads}); do
        part=\$(printf "%03d" \$i)
        touch output/${sra}_1.part_\${part}.fasta.gz
        touch output/${sra}_2.part_\${part}.fasta.gz
    done
    """ : """
    touch output/${sra}_1.fasta
    touch output/${sra}_2.fasta
    """
    """
    mkdir -p output
    ${stub_split_cmd}
    """
}
