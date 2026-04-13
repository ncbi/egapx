#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../utilities'

workflow setup_genome {
    take:
        genome
        organelles
        parameters  // Map : extra parameter and parameter update
    main:
        info = get_genome_info(genome, organelles, true, parameters)
    emit:
        seqid_list = info.seqid_list
        gencoll_asn = info.gencoll_asn
        unpacked_genome = info.fasta
        genome_asn = info.genome_asn
        genome_asnb = info.genome_asnb
        max_intron = info.max_intron
}


workflow setup_ext_genome {
    take:
        genome
        organelles
        parameters  // Map : extra parameter and parameter update
    main:
        info = get_genome_info(genome, organelles, false, parameters)
    emit:
        seqid_list = info.seqid_list
        gencoll_asn = info.gencoll_asn
        unpacked_genome = info.fasta
        genome_asn = info.genome_asn
        genome_asnb = info.genome_asnb
        max_intron = info.max_intron
}


process get_genome_info {
    label 'single_cpu'
    label 'small_mem'
    input:
        path input_genome_file, stageAs: 'src/*'
        path organelles
        val  localize_ids
        val  parameters
    output:
        path '*.seqids', emit: 'seqid_list'
        path '*.asn', emit: 'gencoll_asn'
        path "${out_fasta}", emit: 'fasta'
        path "${genome_asn}", emit: 'genome_asn'
        path "${genome_asnb}", emit: 'genome_asnb'
        env  max_intron, emit: 'max_intron'
    script:
        need_zcat = input_genome_file.toString().endsWith('.gz')
        base_name_stripped = input_genome_file.baseName.toString().replaceAll(/\.(fa(sta)?|fna|asn(t|b)?)(\.gz)?$/, "")
        indexed_fasta_name = input_genome_file.baseName.toString().replaceFirst(/\.(fa(sta)?|fna)(\.gz)?$/, ".fasta")
        if (! indexed_fasta_name.endsWith(".fasta")) {
            indexed_fasta_name += ".fasta"
        }
        genome_dir = "genome"
        fasta_dir = "fasta"
        out_fasta = fasta_dir + "/" + base_name_stripped + ".fasta"
        genome_asn = genome_dir + "/" + base_name_stripped + ".asn"
        genome_asnb = genome_dir + "/" + base_name_stripped + ".asnb"
        max_intron = parameters.max_intron
        genome_size_threshold = parameters.genome_size_threshold
        submitter = parameters.annotation_provider
        name_prefix = parameters.annotation_name_prefix
        lcl_id_cmd = "cat"
        if(localize_ids) {
            lcl_id_cmd = "sed -r 's/>([^ |]+)( .*)?\$/>lcl|\\1\\2/'"
        }
        add_title_cmd = "sed -r 's/>([^ ]+)\$/>\\1 title/'"
    """
    # echo "need_zcat: ${need_zcat}, out_fasta: ${out_fasta}"
    mkdir -p ${genome_dir}
    mkdir -p ${fasta_dir}
    mkdir -p tmp/asncache
    
    first_char=""
    if [[ ${need_zcat} == true ]]; then
        first_char=`zcat $input_genome_file | head -c 100 | grep -oP "[[:print:]]" | head -c 1 || true`
    else
        first_char=`head -c 100 $input_genome_file | grep -oP "[[:print:]]" | head -c 1`
    fi

    if [[ \${first_char} == ">" ]]; then
        if [[ ${need_zcat} == true ]]; then
            zcat ${input_genome_file} | ${lcl_id_cmd} | ${add_title_cmd} > ${out_fasta}
        else
            ${lcl_id_cmd} ${input_genome_file} | ${add_title_cmd} > ${out_fasta}
        fi
        multireader -flags ParseRawID -out-format asn_text -input ${out_fasta} -output ${genome_asn}
        multireader -flags ParseRawID -out-format asn_binary -input ${out_fasta} -output ${genome_asnb}
        prime_cache -cache tmp/asncache/ -ifmt asnb-seq-entry  -i ${genome_asnb} -oseq-ids ./cache_id_list  -split-sequences
    elif [[ \${first_char} == "S" ]]; then
        if [[ ${need_zcat} == true ]]; then
            zcat ${input_genome_file}  > ${genome_asn}
        else
            cp ${input_genome_file} ${genome_asn}
        fi
        asn_translator -i  ${genome_asn} -o ${genome_asnb} -b 
        prime_cache -cache tmp/asncache/ -ifmt asn-seq-entry  -i ${input_genome_file} -oseq-ids ./cache_id_list  -split-sequences
        getfasta -nogenbank -asn-cache tmp/asncache/ -i ./cache_id_list  -o ${out_fasta}
    else 
        if [[ ${need_zcat} == true ]]; then
            zcat ${input_genome_file}  > ${genome_asnb}
        else
            cp ${input_genome_file} ${genome_asnb}
        fi
        asn_translator -i  ${genome_asnb} -o ${genome_asn}  
        prime_cache -cache tmp/asncache/ -ifmt asnb-seq-entry  -i ${input_genome_file} -oseq-ids ./cache_id_list  -split-sequences
        getfasta -nogenbank -asn-cache tmp/asncache/ -i ./cache_id_list  -o ${out_fasta}
    fi

    # Old way, now use gc_get_molecules. For multipart ids with gi first use the second part
    # grep -oP "^>\\K[^ ]+" ${out_fasta} | sed 's/^\\(gi|[0-9]\\+\\)|\\([^|]\\+|[^|]\\+\\)|\\?/\\2/' >list.seqids
    # Using all parts of multipart ids is preferrable, but slower - one more pass over genomic FASTA
    
    ## gc_create seq-id works with some characters that gc_create fasta does not
    gc_create -unplaced ./cache_id_list -unplaced-fmt seq-id -fasta-parse-raw-id -gc-assm-name "$name_prefix" -nogenbank -asn-cache tmp/asncache/ >${base_name_stripped}-gencoll.asn
    
    gc_get_molecules -gc-assembly ${base_name_stripped}-gencoll.asn -filter all -level top-level > list.seqids

    #TODO: subtract organelles from list

    # Get exact genome size using seqkit.
    genome_size=`seqkit stat --tabular ${out_fasta} | tail -n1 | cut -f5`
    # Max intron logic
    if [ $genome_size_threshold -gt 0 ] && [ \$genome_size -lt $genome_size_threshold ]; then
        # scale max intron to genome size, rounding up to nearest 100kb
        (( max_intron = ($max_intron * genome_size / $genome_size_threshold + 99999) / 100000 * 100000 ))
        # echo "Setting max_intron to \$max_intron"
    else
        max_intron=$max_intron
    fi
    # NB: this printout is essential for effective max_intron value to be reported
    echo "max_intron \$max_intron"
    rm -rf tmp
    """
    
    stub:
        base_name_stripped = input_genome_file.baseName.toString().replaceAll(/\.(fa(sta)?|fna)(\.gz)?$/, "")
        indexed_fasta_name = input_genome_file.baseName.toString().replaceFirst(/\.(fa(sta)?|fna)(\.gz)?$/, ".fasta")
        if (! indexed_fasta_name.endsWith(".fasta")) {
            indexed_fasta_name += ".fasta"
        }
        genome_dir = "genome"
        fasta_dir = "fasta"
        out_fasta = fasta_dir + "/" + base_name_stripped + ".fasta"
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
    max_intron=100001
    echo "Processing genome $input_genome_file"
    # NB: this printout is essential for effective max_intron value to be reported
    echo "max_intron \$max_intron"
    """
}


workflow setup_proteins {
    take:
        proteins
        filter_taxons
        additional_proteins
        parameters  // Map : extra parameter and parameter update
    main:
        convert_proteins(proteins, true, filter_taxons, additional_proteins)
    emit:
        unpacked_proteins = convert_proteins.out.unpacked_proteins
        proteins_asn = convert_proteins.out.proteins_asn
        proteins_asnb = convert_proteins.out.proteins_asnb
}


workflow setup_ext_proteins {
    take:
        proteins
        parameters  // Map : extra parameter and parameter update
    main:
        convert_proteins(proteins, false, [], [])
    emit:
        unpacked_proteins = convert_proteins.out.unpacked_proteins
        proteins_asn = convert_proteins.out.proteins_asn
        proteins_asnb = convert_proteins.out.proteins_asnb
}


process convert_proteins {
    label 'single_cpu'
    label 'small_mem'
    input:
        path fasta_proteins_file, stageAs: 'src/*'
        val  localize_ids
        val  filter_taxons
        path additional_proteins, stageAs: 'user_src/*'
    output:
        path out_fasta, emit: 'unpacked_proteins'
        path proteins_asn, emit: 'proteins_asn'
        path proteins_asnb, emit: 'proteins_asnb'
    script:
        def need_zcat = fasta_proteins_file.toString().endsWith('.gz')
        def base_name_stripped = fasta_proteins_file.baseName.toString().replaceAll(/\.(fa(sta)?|faa)(\.gz)?$/, "")
        def fasta_name = base_name_stripped + ".faa"

        asn_dir = "asn"
        fasta_dir = "fasta"
        out_fasta = fasta_dir + "/" + fasta_name
        proteins_asn = asn_dir + "/" + base_name_stripped + ".asn"
        proteins_asnb = asn_dir + "/" + base_name_stripped + ".asnb"
        def source_cmd = need_zcat ? "zcat" : "cat"
        def lcl_id_cmd = localize_ids ? " | sed 's/>\\([^ |]\\+\\)\\( .*\\)\\?\$/>lcl\\|\\1\\2/'" : ""
        def fiter_taxons_cmd = filter_taxons ? ' | seqkit grep -r -n -p "\\[(' + filter_taxons.join(")|(") + ')\\]"' : ""
    """
    mkdir -p ${asn_dir}
    mkdir -p ${fasta_dir}
    mkdir -p tmp/interim
    ${source_cmd} ${fasta_proteins_file}${fiter_taxons_cmd}${lcl_id_cmd} > tmp/interim/out.faa
    for f in user_src/*; do
        [ -e "\$f" ] || continue
        if [[ \$f == *.gz ]]; then
            zcat \$f >> tmp/interim/out.faa
        else
            cat \$f >> tmp/interim/out.faa
        fi
    done
    seqkit rmdup tmp/interim/out.faa -o ${out_fasta}
    multireader -flags ParseRawID -out-format asn_text -input ${out_fasta} -output ${proteins_asn} -parse-mods
    multireader -flags ParseRawID -out-format asn_binary -input ${out_fasta} -output ${proteins_asnb} -parse-mods
    rm -rf tmp
    """

    stub:
        def need_zcat = fasta_proteins_file.toString().endsWith('.gz')
        def base_name_stripped = fasta_proteins_file.baseName.toString().replaceAll(/\.(fa(sta)?|faa)(\.gz)?$/, "")
        def fasta_name = base_name_stripped + ".faa"

        asn_dir = "asn"
        fasta_dir = "fasta"
        out_fasta = fasta_dir + "/" + fasta_name
        proteins_asn = asn_dir + "/" + base_name_stripped + ".asn"
        proteins_asnb = asn_dir + "/" + base_name_stripped + ".asnb"
        def source_cmd = need_zcat ? "zcat" : "cat"
        def lcl_id_cmd = localize_ids ? "| sed 's/>\\([^ |]\\+\\)\\( .*\\)\\?\$/>lcl\\|\\1\\2/'" : ""
        def fiter_taxons_cmd = filter_taxons ? '| seqkit grep -r -n -p "\\[(' + filter_taxons.join(")|(") + ')\\]"' : ""
        print "Protein processing command: ${source_cmd} ${fasta_proteins_file}${fiter_taxons_cmd}${lcl_id_cmd} > ${out_fasta}"
    """
    mkdir -p $asn_dir
    mkdir -p $fasta_dir
    touch $out_fasta
    touch $proteins_asn
    touch $proteins_asnb
    """
}


workflow setup_long_reads {
    take:
        reads
        srr_id
    main:
        rename_fasta_ids(reads, srr_id)
    emit:
        renamed_reads = rename_fasta_ids.out
}


process rename_fasta_ids {
    label 'single_cpu'
    label 'small_mem'
    input:
        tuple val(sampleID), path(fastx, stageAs: "reads/*")
        val  srr_id
    output:
        tuple val(sampleID),  path ('output/*')  , emit: 'fasta_pair_list'
    script:
        def file_name = fastx.getBaseName()
        if (!file_name.endsWith(".fasta")) {
            file_name += ".fasta"
        }
        def srrFmt = String.format('SRR%08d', (srr_id as int))
    """
    mkdir -p output
    seqkit fq2fa -j 7 '${fastx}' | seqkit replace -j 7 -w 0 -p '.*' -r 'gnl|SRA|${srrFmt}.{nr}.1' > output/${file_name}
    """
    stub:
        def file_name = fastx.getBaseName()
        if (!file_name.endsWith(".fasta")) {
            file_name += ".fasta"
        }
    """
    mkdir -p output
    echo $srr_id > output/$file_name
    """
}
