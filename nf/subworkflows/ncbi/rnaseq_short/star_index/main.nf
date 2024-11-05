#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


process build_index {
    debug true
    label 'big_job'
    input:
        path genome_file
        val parameters
    output:
        path "${out_dir}", emit: "index_dir"
    script:
        out_dir = genome_file.toString().replaceFirst(/\.(fa(sta)?|fna)$/, ".index")
    """
    mkdir -p $out_dir
    echo "in $genome_file, out $out_dir"
    # Check that --genomeSAindexNbases is not set
    if [[ ! "$parameters" =~ --genomeSAindexNbases ]]; then
        # Formula from https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf, section 2.2.5
        sa_index_bases=`python3 -c "import os;import math;print(math.floor(min(14, math.log2(os.path.getsize('$genome_file'))/2 - 1)))"`
        effective_parameters="$parameters --genomeSAindexNbases \$sa_index_bases"
    else
        effective_parameters="$parameters"
    fi
    echo "ls"
    ls

    echo "STAR index effective parameters: \$effective_parameters"
    STAR  \$effective_parameters --runMode genomeGenerate --genomeDir $out_dir --genomeFastaFiles $genome_file
    chmod a+rx $out_dir
    STAR  $parameters --runMode genomeGenerate --genomeDir $out_dir --genomeFastaFiles $genome_file
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
        default_params = "--runThreadN 6"
        build_index(genome, merge_params(default_params, parameters, 'STAR'))
   emit:
        index_dir = build_index.out.index_dir
}
