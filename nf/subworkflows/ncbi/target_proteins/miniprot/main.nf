#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


def get_effective_params(parameters, max_intron) {
    def default_params = "-t 8 -G ${max_intron}"
    def value = parameters.get("miniprot", "")
    value = value.replaceFirst("-cpu-count", "-t")
    value = value.replaceFirst("-max-intron", "-G")
    parameters['miniprot'] = value
    def effective_params = merge_params(default_params, parameters, "miniprot")
    return effective_params
}

workflow miniprot {
    take:
        fasta_genome_file  //path: genome fasta file
        fasta_proteins_file  //path: protein fasta file
        max_intron          //int: max intron length
        parameters      // Map : extra parameter and parameter update
    main:
        // println("Miniprot max intron: ${max_intron}")
        def items_per_chunk = merge_params("-n 1000000000", parameters, "split_proteins").replaceFirst("-n ", "").toInteger()
        def protein_chunks
        if (items_per_chunk == 1000000000) {
            protein_chunks = fasta_proteins_file
        } else {
            protein_chunks = split_proteins(fasta_proteins_file, items_per_chunk)
        }
        run_miniprot(fasta_genome_file, protein_chunks.flatten(), max_intron, parameters)

    emit:
        miniprot_file = run_miniprot.out.miniprot_file
}


process split_proteins {
    input:
        path fasta_proteins_file, stageAs: 'inputs/input_prots.faa*'
        val  items_per_chunk
    output:
        path 'output/*'
    script:
    """
    #!/usr/bin/env python3
    import os

    os.makedirs("output", exist_ok=True)
    with open("${fasta_proteins_file}", 'rt') as f:
        items = 0
        chunk = []
        nextfile = 1
        for line in f:
            if line and line[0] == '>':
                items += 1
                if items >= ${items_per_chunk}:
                    with open(f"output/{nextfile}.prots.faa", "w") as outf:
                        outf.write(''.join(chunk))
                        chunk = []
                        nextfile += 1
                        items = 1
            chunk.append(line)
        if chunk:
            with open(f"output/{nextfile}.prots.faa", "w") as outf:
                outf.write(''.join(chunk))
    """
    stub:
    print("items_per_chunk ${items_per_chunk}")
    """
    mkdir -p output
    touch output/1.fa
    touch output/2.fa
    touch output/3.fa
    """
}


process run_miniprot {
    label 'huge_job'
    label 'long_job'
    input:
        path fasta_genome_file
        path fasta_proteins_file
        val max_intron
        val parameters
    output:
        path ('output/*.paf'), emit: 'miniprot_file'

    script:
        def paf_name = fasta_proteins_file.baseName.toString() + ".paf"
        def effective_params = get_effective_params(parameters, max_intron)
        // println("Miniprot params: ${effective_params}")
    """
    mkdir -p output
    miniprot ${effective_params}  ${fasta_genome_file} ${fasta_proteins_file} > output/${paf_name}
    """
    stub:
        def paf_name = fasta_proteins_file.baseName.toString() + ".paf"
        def effective_params = get_effective_params(parameters, max_intron)
        println("Miniprot params: ${effective_params}")
    """
    mkdir -p output
    touch output/${paf_name}
    """
}
