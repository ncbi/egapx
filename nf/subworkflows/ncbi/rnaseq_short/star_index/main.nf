#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


process build_index {
    label 'big_job'
    input:
        path genome_file
        val parameters
    output:
        path out_dir
    script:
        out_dir = genome_file.toString().replaceFirst(/\.(fa(sta)?|fna)$/, ".index")
    """
    echo "in $genome_file, out $out_dir"
    STAR  $parameters --runMode genomeGenerate --genomeDir $out_dir --genomeFastaFiles $genome_file
    chmod a+rx $out_dir
    """

    stub:
        out_dir = genome_file.toString().replaceFirst(/\.(fa(sta)?|fna)$/, ".index")
    """
    mkdir -p ${out_dir}
    touch ${out_dir}/Genome
    touch ${out_dir}/Log.out
    touch ${out_dir}/SA
    touch ${out_dir}/SAindex
    touch ${out_dir}/chrLength.txt
    touch ${out_dir}/chrName.txt
    touch ${out_dir}/chrNameLength.txt
    touch ${out_dir}/chrStart.txt
    touch ${out_dir}/genomeParameters.txt
    """
}


workflow star_index {
   take:
        genome      //path: genome
        parameters  // Map : extra parameter and parameter update
   main:
        default_params = "--runThreadN 6 --genomeSAindexNbases 14"
        build_index(genome, merge_params(default_params, parameters, 'STAR'))
   emit:
        build_index.out
}
