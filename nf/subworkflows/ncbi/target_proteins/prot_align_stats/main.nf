#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow prot_align_stats {
    take:
        align_tab
        assembly_taxid_list
        proteins_accessions_list
        filter_taxons
        annotation_name_prefix
        parameters  // Map : extra parameter and parameter update
    main:
        String stats_params =  merge_params("", parameters, 'prot_align_stats')
        run_prot_align_stats(align_tab, assembly_taxid_list, proteins_accessions_list, filter_taxons, annotation_name_prefix, stats_params)
    emit:
        prot_align_stats = run_prot_align_stats.out.prot_align_stats
}


/*
process run_prot_align_stats {
    label 'single_cpu'
    label 'small_mem'
    input:
        path gencoll_asn
        path proteins_asn
        path prosplign_aligns
        path aligns_to_gnomon
        path bags
        val taxid
        val parameters
    output:
        path "prot_align_stats.xml", emit: "prot_align_stats"
    script:
    """
    mkdir -p tmp/asncache
    auto_prime_cache.py -cache tmp/asncache/ -i ${proteins_asn} -oseq-ids spids2 -split-sequences
    echo "${aligns_to_gnomon.join('\n')}" > ./aligns_to_gnomon.mft
    echo "${prosplign_aligns.join('\n')}" > ./prosplign_aligns.mft        
    echo "${bags.join('\n')}" > ./bag_sources.mft
    /netmnt/vast01/gpi/regr/GPIPE_REGR1/system/current/arch/x86_64/bin/protein_align_stats -nogenbank -asn-cache tmp/asncache/ -gencoll-asn $gencoll_asn -skip-pig-lookup $parameters -bags bag_sources.mft  -prosplign-aligns prosplign_aligns.mft -aligns-to-gnomon aligns_to_gnomon.mft -taxid $taxid -output prot_align_stats.xml
    """
    stub:
    """    touch prot_align_stats.xml
    """
}

*/

process run_prot_align_stats {
    label 'single_cpu'
    label 'small_mem'
    input:
        path align_tab
        path assembly_taxid_list
        path proteins_accessions_list
        val filter_taxons
        val annotation_name_prefix
        val parameters
    output:
        path "./prot_align_stats.xml", emit: "prot_align_stats"
    script:
        def filter_taxons_str = "[\"${filter_taxons.join('","')}\"]"
    """
    #!/usr/bin/env python3

    from pathlib import Path
    import xml.etree.ElementTree as ET

    assembly_taxid_list = '${assembly_taxid_list}'
    proteins_accessions_list = '${proteins_accessions_list}'
    align_tab = '${align_tab}'

    data_dict2 = {}

    reverse_data_dict = {}
    taxa = ${filter_taxons_str}
    try:
        with open(assembly_taxid_list, mode='r', encoding='utf-8') as f:
            for line in f:
                # Split by tab character
                parts = line.strip().split('\\t')
                if parts:
                    # key is first column, value is the rest of the data
                    taxid_list = parts[2:][0].strip().split(';')  # Get all columns after the first two
                    data_dict2[parts[1]] = parts[0]
                    reverse_data_dict[parts[0]] = parts[1]
    except Exception as e:
        print(f"Error reading assembly_taxid_list: {e}")
        open("prot_align_stats.xml", "w").close()  # create empty output file in case of error
        exit(0)

    translation_to_taxid = {}
    taxid_list = []
    if taxa:
        for taxon in taxa:
            translation_to_taxid[taxon] = data_dict2.get(taxon, [])
            taxid_list.append(data_dict2.get(taxon, []))
            print(f"Taxa from run parameters: {taxa}") 

    accessions_taxid = {}
    try:
         with open(proteins_accessions_list, mode='r', encoding='utf-8') as file:
            next(file)  # Skip the header line
            for line in file:
                a = line.strip().split('\\t')
                accessions = a[1]
                taxid = a[2]
                if taxid in taxid_list:
                    #print(f"Taxid {taxid} is in the unique taxid list, adding accession {accessions}")
                    res = accessions_taxid.get(accessions , None)
                    if res is not None:
                        print(f"Accession {accessions} already exists in the dictionary with taxid {res}, skipping, current taxid is {taxid}")
                    else:
                        #print(f"Adding accession {accessions} with taxid {taxid} to the dictionary")
                        accessions_taxid[accessions] = taxid
    except Exception as e:
        print(f"Error reading proteins_accessions_list: {e}")
        open("prot_align_stats.xml", "w").close()  # create empty output file in case of error
        exit(0)
    
    #print(accessions_taxid)
    acc_taxid = {}
    for item in accessions_taxid.items():
    #   print(f"Accession: {item[0]}, Taxid: {item[1]}")
        acc = item[0][0:3]
        taxid = item[1]
        ret = acc_taxid.get(acc, None)
        if ret is None:
            acc_taxid[acc] = {taxid : 1}
        else:
            acc_taxid[acc][taxid] = acc_taxid[acc].get(taxid, 0) + 1

    #print(acc_taxid)

    acc_taxid_ident_cov = {}
    try:
         with open(align_tab, mode='r', encoding='utf-8') as file:
            next(file)  # Skip the header line
            for line in file:
                parts = line.strip().split('\\t')
                if parts:
                    acc = parts[0][0:3]
                    taxid = parts[1]
                    ident = float(parts[2])
                    coverage = float(parts[3])
                    #print(f"Accession: {acc}, Taxid: {taxid}, Identity: {ident}, Coverage: {coverage}")
                    ret = acc_taxid_ident_cov.get(acc, None)
                    if ret is None:
                        acc_taxid_ident_cov[acc] = {taxid: {'ident': ident, 'coverage': coverage, 'count': 1}}
                    else:
                        if taxid not in acc_taxid_ident_cov[acc]:
                            acc_taxid_ident_cov[acc][taxid] = {'ident': ident, 'coverage': coverage, 'count': 1}
                        else:
                            acc_taxid_ident_cov[acc][taxid]['count'] += 1
                            acc_taxid_ident_cov[acc][taxid]['ident'] += ident
                            acc_taxid_ident_cov[acc][taxid]['coverage'] += coverage
    except Exception as e:
        print(f"Error reading align_tab: {e}")
        open("prot_align_stats.xml", "w").close()  # create empty output file in case of error
        exit(0)

    #print(acc_taxid_ident_cov)

    combined_data = {}
    for acc, taxid_info in acc_taxid_ident_cov.items():
        for taxid, ident_cov in taxid_info.items():
            ret = acc_taxid.get(acc, None)
            total_count = 0
            if ret is not None:
                total_count = acc_taxid[acc].get(taxid, 0)
                if total_count != 0:
                    del acc_taxid[acc][taxid]  # Remove the count from the previous dictionary to avoid double counting
            align_count = ident_cov['count']
            total_count += align_count
            avg_ident = ident_cov['ident'] / align_count
            avg_cov = ident_cov['coverage'] / align_count

            combined_data[(acc, taxid)] = { 'percent aligned': 100 * align_count / total_count, 'avg_ident': avg_ident, 'avg_cov': avg_cov, 'species name': reverse_data_dict.get(taxid, 'Unknown'), 'AlignCount': align_count, 'TotalSequences': total_count}

    for acc, taxid_count in acc_taxid.items():
        for taxid, count in taxid_count.items():
            if (acc, taxid) not in combined_data:
                combined_data[(acc, taxid)] = {'percent aligned': 0, 'avg_ident': 0, 'avg_cov': 0, 'species name': reverse_data_dict.get(taxid, 'Unknown'), 'AlignCount': 0, 'TotalSequences': count}  



    root = ET.Element("ProteinAlignStats")
    assAcc = ET.SubElement(root, "AssemblyStats", accession="$annotation_name_prefix")
        
    for key, value in combined_data.items():
        acc_prefix = key[0]
        taxid = key[1]
        percent_aligned = value['percent aligned']
        avg_ident = value['avg_ident']
        avg_cov = value['avg_cov']
        species_name = value['species name']
        total_sequences = value['TotalSequences']
        align_count = value['AlignCount']
        category_stats = ET.SubElement(assAcc, "CategoryStats", Species=species_name, Taxid=taxid, Category=acc_prefix)
        ET.SubElement(category_stats, "TotalSequences").text = str(total_sequences)
        ET.SubElement(category_stats, "Aligned").text = str(align_count)
        ET.SubElement(category_stats, "PassedToGnomon").text = str(align_count)
        ET.SubElement(category_stats, "PercentAligned").text = str(int(round(percent_aligned, 0))) + "%"
        ET.SubElement(category_stats, "MeanPctIdentity").text = str(int(round(avg_ident, 0))) + "%"
        ET.SubElement(category_stats, "MeanPctCoverage").text = str(int(round(avg_cov, 0))) + "%"

    ET.indent(root, space="  ", level=0)

    tree = ET.ElementTree(root)
    tree.write("prot_align_stats.xml", encoding="utf-8", xml_declaration=True)

    """
    stub:
    """    touch prot_align_stats.xml
    """
}
