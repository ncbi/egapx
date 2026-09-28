#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'
include { align_as_tabular_sa as align_as_tabular_sa_best } from '../align_as_tabular_sa/main'
include { align_as_tabular_sa as align_as_tabular_sa_filtered } from '../align_as_tabular_sa/main'


workflow prot_align_stats {
    take:
        genome_asn                      // path: genome asn file, needed to tabulate best_aligned_prot_asn/filtered_protein_alignments_asn
        proteins_asn                    // path: proteins asn file, needed to tabulate best_aligned_prot_asn/filtered_protein_alignments_asn
        align_tab
        assembly_taxid_list
        proteins_accessions_list
        filter_taxons
        annotation_name_prefix
        best_aligned_prot_asn            // path: best_aligned_prot alignments, raw asn
        filtered_protein_alignments_asn  // path: filtered alignments, raw asn
        aligner_name
        tabular_conversion_params        // Map : extra parameter and parameter update for the asn-to-tabular conversion
        parameters  // Map : extra parameter and parameter update
    main:
        align_as_tabular_sa_best(genome_asn, proteins_asn, best_aligned_prot_asn, tabular_conversion_params)
        def best_aligned_prot_tab = align_as_tabular_sa_best.out.align_tab
        align_as_tabular_sa_filtered(genome_asn, proteins_asn, filtered_protein_alignments_asn, tabular_conversion_params)
        def filtered_protein_alignments_tab = align_as_tabular_sa_filtered.out.align_tab

        String stats_params =  merge_params("", parameters, 'prot_align_stats')
        run_prot_align_stats(align_tab, assembly_taxid_list, proteins_accessions_list, filter_taxons, annotation_name_prefix, best_aligned_prot_tab, filtered_protein_alignments_tab, aligner_name, stats_params)
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
        path best_aligned_prot_tab, stageAs: 'best_aligned_prot/align.tab'
        path filtered_protein_alignments_tab, stageAs: 'filtered_protein_alignments/align.tab'
        val aligner_name
        val parameters
    output:
        path "./prot_align_stats.xml", emit: "prot_align_stats"
    script:
        def filter_taxons_str = "[\"${filter_taxons.join('","')}\"]"
    """
    #!/usr/bin/env python3

    from collections import Counter
    import re
    import xml.etree.ElementTree as ET

    assembly_taxid_list = '${assembly_taxid_list}'
    proteins_accessions_list = '${proteins_accessions_list}'
    align_tab = '${align_tab}'
    best_aligned_prot_tab = '${best_aligned_prot_tab}'
    filtered_protein_alignments_tab = '${filtered_protein_alignments_tab}'

    aligned_tag = {'miniprot': 'AlignedByMiniprot', 'prosplign': 'AlignedByProSplign'}.get('${aligner_name}', 'Aligned')

    # data_dict2: taxon name -> taxid, reverse_data_dict: taxid -> taxon name
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
                    data_dict2[parts[1]] = parts[0]
                    reverse_data_dict[parts[0]] = parts[1]
    except Exception as e:
        print(f"Error reading assembly_taxid_list: {e}")
        open("prot_align_stats.xml", "w").close()  # create empty output file in case of error
        exit(0)

    # Translate the requested taxa names into taxids; stats are reported only for these
    taxid_list = []
    if taxa:
        for taxon in taxa:
            taxid_list.append(data_dict2.get(taxon, []))
            print(f"Taxa from run parameters: {taxa}") 

    version_re = re.compile(r'[0-9]*[.][0-9]*')

    def strip_version(accession):
        # Mirrors: sed 's/[0-9]*[.][0-9]*//' -- drops the version suffix, leaving the accession prefix
        return version_re.sub('', accession, count=1)

    def count_total_sequences(path):
        # Mirrors: cut -f2-3 | grep '[_]' | sed 's/[0-9]*[.][0-9]*//' | sort | uniq -c
        counts = Counter()
        try:
            with open(path, mode='r', encoding='utf-8') as f:
                next(f)  # Skip the header line
                for line in f:
                    parts = line.rstrip('\\n').split('\\t')
                    if len(parts) < 3:
                        continue
                    accession, taxid = parts[1], parts[2]
                    if '_' not in accession:
                        continue
                    counts[(strip_version(accession), taxid)] += 1
        except Exception as e:
            print(f"Error reading proteins_accessions_list: {e}")
            open("prot_align_stats.xml", "w").close()  # create empty output file in case of error
            exit(0)
        return counts

    def count_unique_pairs(path):
        # Mirrors: grep -v '^#' | cut -f1-2 | sort | uniq | grep '[_]' | sed 's/[0-9]*[.][0-9]*//' | sort | uniq -c
        pairs = set()
        try:
            with open(path, mode='r', encoding='utf-8') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.rstrip('\\n').split('\\t')
                    if len(parts) < 2:
                        continue
                    pairs.add((parts[0], parts[1]))
        except Exception as e:
            print(f"Error reading {path}: {e}")
            open("prot_align_stats.xml", "w").close()  # create empty output file in case of error
            exit(0)
        counts = Counter()
        for accession, taxid in pairs:
            if '_' not in accession:
                continue
            counts[(strip_version(accession), taxid)] += 1
        return counts

    total_sequences_counts = count_total_sequences(proteins_accessions_list)
    aligned_counts = count_unique_pairs(best_aligned_prot_tab)
    passed_to_gnomon_counts = count_unique_pairs(filtered_protein_alignments_tab)

    # acc_taxid_ident_cov: accumulates identity/coverage sums and counts per (prefix, taxid), for later averaging
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

    # Merge the three counters per (prefix, taxid), keeping only the requested taxa
    combined_data = {}
    prefix_taxid_pairs = set(total_sequences_counts) | set(aligned_counts) | set(passed_to_gnomon_counts)
    for acc, taxid in prefix_taxid_pairs:
        if taxid not in taxid_list:
            continue
        total_sequences = total_sequences_counts.get((acc, taxid), 0)
        align_count = aligned_counts.get((acc, taxid), 0)
        passed_to_gnomon = passed_to_gnomon_counts.get((acc, taxid), 0)
        percent_aligned = 100 * align_count / total_sequences if total_sequences else 0
        percent_filter = 100 * passed_to_gnomon / total_sequences if total_sequences else 0
        ident_cov = acc_taxid_ident_cov.get(acc, {}).get(taxid, {})
        ident_cov_count = ident_cov.get('count', 0)
        avg_ident = ident_cov['ident'] / ident_cov_count if ident_cov_count else 0
        avg_cov = ident_cov['coverage'] / ident_cov_count if ident_cov_count else 0

        combined_data[(acc, taxid)] = {'percent aligned': percent_aligned, 'percent filter': percent_filter, 'avg_ident': avg_ident, 'avg_cov': avg_cov, 'species name': reverse_data_dict.get(taxid, 'Unknown'), 'AlignCount': align_count, 'PassedToGnomon': passed_to_gnomon, 'TotalSequences': total_sequences}

    # Emit one CategoryStats element per (accession prefix, taxid) combination
    root = ET.Element("ProteinAlignStats")
    assAcc = ET.SubElement(root, "AssemblyStats", accession="$annotation_name_prefix")
        
    for key, value in combined_data.items():
        acc_prefix = key[0]
        taxid = key[1]
        percent_aligned = value['percent aligned']
        percent_filter = value['percent filter']
        avg_ident = value['avg_ident']
        avg_cov = value['avg_cov']
        species_name = value['species name']
        total_sequences = value['TotalSequences']
        align_count = value['AlignCount']
        passed_to_gnomon = value['PassedToGnomon']
        category_stats = ET.SubElement(assAcc, "CategoryStats", Species=species_name, Taxid=taxid, Category=acc_prefix)
        ET.SubElement(category_stats, "TotalSequences").text = str(total_sequences)
        ET.SubElement(category_stats, aligned_tag).text = str(align_count)
        ET.SubElement(category_stats, "PassedFilter").text = str(passed_to_gnomon)
        ET.SubElement(category_stats, "PercentFilter").text = str(int(percent_filter + 0.5)) + "%"  # round half up, not banker's rounding
        ET.SubElement(category_stats, "PercentAligned").text = str(int(percent_aligned + 0.5)) + "%"  # round half up, not banker's rounding
        ET.SubElement(category_stats, "MeanPctIdentity").text = str(int(avg_ident + 0.5)) + "%"  # round half up, not banker's rounding
        ET.SubElement(category_stats, "MeanPctCoverage").text = str(int(avg_cov + 0.5)) + "%"  # round half up, not banker's rounding

    ET.indent(root, space="  ", level=0)

    tree = ET.ElementTree(root)
    tree.write("prot_align_stats.xml", encoding="utf-8", xml_declaration=True)

    """
    stub:
    """    touch prot_align_stats.xml
    """
}
