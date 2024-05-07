#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/*
 *Execution of: 
 * /netmnt/vast01/gpi/prod/GPIPE_PROD/system/current/arch/x86_64/bin/gp_makeblastdb 
    -nogenbank
    -asn-cache ''
    -db ./spdb
    -asnb swissprot.asnb
    -dbtype prot
    -title 'BLASTdb used for naming tasks in GPipe'

 * /netmnt/vast01/gpi/regr/GPIPE_REGR1/system/2024-03-27.prod.build25780/bin/align_sort 
 *   -ifmt seq-align-set -input-manifest /netmnt/vast01/gpi/regr/GPIPE_REGR1/data00/Gavia_stellata/GP37025.85624/846757/protein_filter.8201942/inp/blast_align.mft 
 *   -k query,query_start,-query_end,query_strand,subject,subject_start,-subject_end,subject_strand,-num_ident,gap_count 
 *   -o /netmnt/vast01/gpi/regr/GPIPE_REGR1/data00/Gavia_stellata/GP37025.85624/846757/protein_filter.8201942/var/blast_align_sorted.asnb 
 *   -tmp /netmnt/vast01/gpi/regr/GPIPE_REGR1/data00/Gavia_stellata/GP37025.85624/846757/protein_filter.8201942/var 
 *   -nogenbank -limit-mem 13G 
 * /netmnt/vast01/gpi/regr/GPIPE_REGR1/system/2024-03-27.prod.build25780/bin/protein_filter 
 *   -lds2 /netmnt/vast01/gpi/regr/GPIPE_REGR1/data00/Gavia_stellata/GP37025.85624/846757/prot_gnomon_prepare.8202002/out/LDS2 
 *   -asn-cache /netmnt/vast01/gpi/regr/GPIPE_REGR1/data00/Gavia_stellata/GP37025.85624/sequence_cache 
 *   -filter 'pct_coverage >= 50' 
 *   -input /netmnt/vast01/gpi/regr/GPIPE_REGR1/data00/Gavia_stellata/GP37025.85624/846757/protein_filter.8201942/var/blast_align_sorted.asnb 
 *   -nr-path /netmnt/vast01/gpi/regr/GPIPE_REGR1/data00/Gavia_stellata/GP37025.85624/846757/create_nr_blast_db.8200792/out/blastdb 
 *   -output /netmnt/vast01/gpi/regr/GPIPE_REGR1/data00/Gavia_stellata/GP37025.85624/846757/protein_filter.8201942/out/blast.asn 
 *   -taxid 37040 
 */

include { merge_params; to_map; shellSplit } from '../../utilities'


workflow best_protein_hits {
    take:
        gnomon_prot_asn
        swiss_prot_asn
        prot_alignments
        taxid
        parameters  // Map : extra parameter and parameter update
    main:
        String align_sort_params = merge_params('-nogenbank', parameters, 'align_sort')
        String protein_filter_params = merge_params('-nogenbank', parameters, 'protein_filter')

        run_protein_filter(gnomon_prot_asn, swiss_prot_asn, prot_alignments, taxid, align_sort_params, protein_filter_params)

    emit:
        alignments = run_protein_filter.out
}



process run_protein_filter {
    input:
        path gnomon_prot_asn, stageAs: 'indexed/*'
        path swiss_prot_asn, stageAs: 'indexed/*'
        path input_prot_alignments, stageAs: "input_alignments.asnb"
        val taxid
        val align_sort_params
        val protein_filter_params
    output:
        path "output/*"
    script:
    """
    
    mkdir -p ./asncache/

    prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ${gnomon_prot_asn} -oseq-ids /dev/null -split-sequences
    prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ${swiss_prot_asn} -oseq-ids /dev/null -split-sequences

    mkdir ./output
    mkdir ./work
    mkdir ./spdb

    /netmnt/vast01/gpi/prod/GPIPE_PROD/system/current/arch/x86_64/bin/gp_makeblastdb -nogenbank -asn-cache ./asncache/ -db ./spdb/spdb -asnb ${swiss_prot_asn} -dbtype prot -title 'BLASTdb used for naming tasks in GPipe'

    align_sort $align_sort_params -asn-cache ./asncache/ -input $input_prot_alignments -o ./sorted_prot_alignments.asnb  -tmp ./work/ 

    protein_filter $protein_filter_params -asn-cache ./asncache/ -taxid $taxid -input ./sorted_prot_alignments.asnb -output ./output/best_protein_hits.asnb -nr-path ./spdb/spdb
    
    rm -rf ./work
    """

    stub:
    """
    mkdir -p output
    touch ./output/best_protein_hits.asnb
    """
}


