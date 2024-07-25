#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../utilities'

workflow setup_genome {
    take:
        genome
        organelles
        parameters  // Map : extra parameter and parameter update
    main:
        get_genome_info(genome, organelles, parameters)
    emit:
        seqid_list = get_genome_info.out.seqid_list
        gencoll_asn = get_genome_info.out.gencoll_asn
        unpacked_genome = get_genome_info.out.fasta
        genome_asn = get_genome_info.out.genome_asn
        genome_asnb = get_genome_info.out.genome_asnb
        max_intron = get_genome_info.out.max_intron
}


process get_genome_info {
    debug true
    input:
        path fasta_genome_file, stageAs: 'src/*'
        path organelles
        val  parameters
    output:
        path '*.seqids', emit: 'seqid_list'
        path '*.asn', emit: 'gencoll_asn'
        path "${out_fasta}", emit: 'fasta'
        path "${genome_asn}", emit: 'genome_asn'
        path "${genome_asnb}", emit: 'genome_asnb'
        env  max_intron, emit: 'max_intron'
    script:
        need_zcat = fasta_genome_file.toString().endsWith('.gz')
        base_name_stripped = fasta_genome_file.baseName.toString().replaceAll(/\.(fa(sta)?|fna)(\.gz)?$/, "")
        indexed_fasta_name = fasta_genome_file.baseName.toString().replaceFirst(/\.(fa(sta)?|fna)(\.gz)?$/, ".fasta")
        if (! indexed_fasta_name.endsWith(".fasta")) {
            indexed_fasta_name += ".fasta"
        }
        genome_dir = "genome"
        fasta_dir = "fasta"
        out_fasta = fasta_dir + "/" + indexed_fasta_name
        genome_asn = genome_dir + "/" + base_name_stripped + ".asn"
        genome_asnb = genome_dir + "/" + base_name_stripped + ".asnb"
        max_intron = parameters.max_intron
        genome_size_threshold = parameters.genome_size_threshold        
    """
    # echo "need_zcat: ${need_zcat}, out_fasta: ${out_fasta}"
    mkdir -p ${genome_dir}
    mkdir -p ${fasta_dir}
    if [[ ${need_zcat} == true ]]; then
        # zcat ${fasta_genome_file} | sed 's/^\\(>gi|[0-9]\\+\\)|\\?\\([^ ]\\+\\)\\(.*\\)/\\1\\3/' > ${out_fasta}
        # zcat ${fasta_genome_file} > ${out_fasta}
        zcat ${fasta_genome_file} | sed 's/>\\([^ |]\\+\\)\\( .*\\)\\?\$/>lcl\\|\\1\\2/' > ${out_fasta}
    else
        # sed 's/^\\(>gi|[0-9]\\+\\)|\\?\\([^ ]\\+\\)\\(.*\\)/\\1\\3/' ${fasta_genome_file} > ${out_fasta}
        # mv ${fasta_genome_file} ${out_fasta}
        sed 's/>\\([^ |]\\+\\)\\( .*\\)\\?\$/>lcl\\|\\1\\2/' ${fasta_genome_file} > ${out_fasta}
    fi
    # Old way, now use gc_get_molecules. For multipart ids with gi first use the second part
    # grep -oP "^>\\K[^ ]+" ${out_fasta} | sed 's/^\\(gi|[0-9]\\+\\)|\\([^|]\\+|[^|]\\+\\)|\\?/\\2/' >list.seqids
    multireader -flags ParseRawID -out-format asn_text -input ${out_fasta} -output ${genome_asn}
    multireader -flags ParseRawID -out-format asn_binary -input ${out_fasta} -output ${genome_asnb}
    lds2_indexer -source ${genome_dir}/ -db LDS2
    # Using all parts of multipart ids is preferrable, but slower - one more pass over genomic FASTA
    gc_create -unplaced ${out_fasta} -unplaced-fmt fasta -fasta-parse-raw-id -gc-assm-name "EGAPx Test Assembly" -nogenbank -lds2 LDS2 >gencoll.asn
    gc_get_molecules -gc-assembly gencoll.asn -filter all -level top-level > list.seqids

    #TODO: subtract organelles from list

    # This is a rough estimate because we don't need the more accurate size
    genome_size=`wc -c <${out_fasta}`
    # Max intron logic
    if [ $genome_size_threshold -gt 0 ] && [ \$genome_size -lt $genome_size_threshold ]; then
        # scale max intron to genome size, rounding up to nearest 100kb
        (( max_intron = ($max_intron * genome_size / $genome_size_threshold + 99999) / 100000 * 100000 ))
        # echo "Setting max_intron to \$max_intron"
    else
        max_intron=$max_intron
    fi
    """
    
    stub:
        base_name_stripped = fasta_genome_file.baseName.toString().replaceAll(/\.(fa(sta)?|fna)(\.gz)?$/, "")
        indexed_fasta_name = fasta_genome_file.baseName.toString().replaceFirst(/\.(fa(sta)?|fna)(\.gz)?$/, ".fasta")
        if (! indexed_fasta_name.endsWith(".fasta")) {
            indexed_fasta_name += ".fasta"
        }
        genome_dir = "genome"
        fasta_dir = "fasta"
        out_fasta = fasta_dir + "/" + indexed_fasta_name
        genome_asn = genome_dir + "/" + base_name_stripped + ".asn"
        genome_asnb = genome_dir + "/" + base_name_stripped + ".asnb"
    """
    mkdir -p $genome_dir
    mkdir -p $fasta_dir
    touch $out_fasta
    touch $genome_asn
    touch $genome_asnb
    touch gencoll.asn
    touch list.seqids
    max_intron=10000
    echo "Processing genome $fasta_genome_file"
    echo "Setting max_intron to \$max_intron"
    """
}


workflow setup_proteins {
    take:
        proteins
        parameters  // Map : extra parameter and parameter update
    main:
        convert_proteins(proteins)
    emit:
        unpacked_proteins = convert_proteins.out.unpacked_proteins
        proteins_asn = convert_proteins.out.proteins_asn
        proteins_asnb = convert_proteins.out.proteins_asnb
}


process convert_proteins {
    input:
        path fasta_proteins_file, stageAs: 'src/*'
    output:
        path out_fasta, emit: 'unpacked_proteins'
        path proteins_asn, emit: 'proteins_asn'
        path proteins_asnb, emit: 'proteins_asnb'
    script:
        need_zcat = fasta_proteins_file.toString().endsWith('.gz')
        base_name_stripped = fasta_proteins_file.baseName.toString().replaceAll(/\.(fa(sta)?|faa)(\.gz)?$/, "")
        fasta_name = base_name_stripped + ".faa"

        asn_dir = "asn"
        fasta_dir = "fasta"
        out_fasta = fasta_dir + "/" + fasta_name
        proteins_asn = asn_dir + "/" + base_name_stripped + ".asn"
        proteins_asnb = asn_dir + "/" + base_name_stripped + ".asnb"
    """
    mkdir -p ${asn_dir}
    mkdir -p ${fasta_dir}
    if [[ ${need_zcat} == true ]]; then
        zcat ${fasta_proteins_file} | sed 's/>\\([^ |]\\+\\)\\( .*\\)\\?\$/>lcl\\|\\1\\2/' > ${out_fasta}
    else
        sed 's/>\\([^ |]\\+\\)\\( .*\\)\\?\$/>lcl\\|\\1\\2/' ${fasta_proteins_file} > ${out_fasta}
    fi
    multireader -flags ParseRawID -out-format asn_text -input ${out_fasta} -output ${proteins_asn}
    multireader -flags ParseRawID -out-format asn_binary -input ${out_fasta} -output ${proteins_asnb}
    """

    stub:
        base_name_stripped = fasta_proteins_file.baseName.toString().replaceAll(/\.(fa(sta)?|faa)(\.gz)?$/, "")
        fasta_name = base_name_stripped + ".faa"
        asn_dir = "asn"
        fasta_dir = "fasta"
        out_fasta = fasta_dir + "/" + fasta_name
        proteins_asn = asn_dir + "/" + base_name_stripped + ".asn"
        proteins_asnb = asn_dir + "/" + base_name_stripped + ".asnb"
    
    """
    mkdir -p $asn_dir
    mkdir -p $fasta_dir
    touch $out_fasta
    touch $proteins_asn
    touch $proteins_asnb
    """

}
