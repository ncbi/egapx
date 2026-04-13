#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/*
 *Execution of:
 * gbin/diamond 
 *    -asn-cache GP37025.85624/sequence_cache 
 *    -blastp-args '--sam-query-len --comp-based-stats 0 --evalue 0.0001 --very-sensitive --max-hsps 3' 
     -diamond-executable gthird-party/diamond/diamond 
 *    -lds2 prot_gnomon_prepare.8202002/out/LDS2 
 *    -ofmt seq-align-set 
 *    -output-dir diamond.8202022/out 
 *    -output-manifest diamond.8202022/out/align.mft 
 *    -output-prefix hits 
 *    ## query is gnomon-made proteins 'gnl|GNOMON|23016146.p'
 *    ## query-fmt is <String, `fasta', `seq-ids'>
 *    -query-fmt seq-ids 
 *    -query-manifest diamond.8202022/inp/query_ids.mft 
 *    ## subject is swiss-prot ids 'sp|A0A009IHW8.1|ABTIR_ACIB9'
 *    -subject-fmt seq-ids 
 *    -subject-manifest diamond.8202022/inp/subject_ids.mft 
 *    -work-area diamond.8202022/tmp

 */

include {to_map; shellSplit } from '../../utilities'


process fetch_swiss_prot_asn  {
    label 'single_cpu'
    label 'small_mem'
    input:
        path swiss_prot_url
    output:
        path "swissprot.asnb", emit: "swiss_prot_asn"
    script:
    """
        zcat ${swiss_prot_url} > swissprot.asnb
    """   
    stub:
    """
        touch swissprot.asnb
    """
}

process get_swiss_prot_ids {
    label 'single_cpu'
    label 'small_mem'
    label 'sqlite3'
    input:
        path swiss_prot_asn, stageAs: 'input/*'
    output:
        path "output/swiss_prot_ids"
    script: 
    """
        mkdir -p output
        mkdir -p tmp
        lds2_indexer  -db tmp/lds -source input
        sqlite3 tmp/lds "SELECT txt_id FROM seq_id WHERE orig=1 AND int_id IS NULL;" > output/swiss_prot_ids
        rm -rf tmp
    """
    stub:
    """
        mkdir -p output
        touch output/swiss_prot_ids
    """
}

process run_diamond_egap {
    label 'multi_cpu'
    label 'med_mem'
    input:
        path gnomon_prot_ids
        path swiss_prot_ids
        path gnomon_prot_asn, stageAs: 'indexed/*'
        path swiss_prot_asn, stageAs: 'indexed/*'
        val params
    output:
        path "output/*"
    script:
    // print(params)
    """
   
    ###diamond_bin=`which diamond`
    #diamond_egap uses GP_HOME to build paths to both some gp apps, and third-party
    #GP_HOME needs to be the directory that contains third-party, and the directory that contains bin/<gp apps> 
    diamond_bin=\${GP_HOME}/third-party/diamond/diamond
    
    mkdir -p tmp/asncache
    if [[ -n \${TMPDIR-} ]]; then
        mkdir -p \${TMPDIR} || true
    fi
    auto_prime_cache.py -cache tmp/asncache/ -i ${gnomon_prot_asn} -oseq-ids /dev/null -split-sequences
    auto_prime_cache.py -cache tmp/asncache/ -i ${swiss_prot_asn} -oseq-ids /dev/null -split-sequences

    mkdir -p output
    mkdir -p tmp/work

    echo  ${params}
    echo "${gnomon_prot_ids.join('\n')}" > query.mft
    diamond_egap  ${params}  -asn-cache tmp/asncache/ -nogenbank -query-manifest query.mft -subject ${swiss_prot_ids} -output-dir ./output/ -work-area tmp/work/  -diamond-executable \${diamond_bin}
    rm -rf tmp
    """

    stub:
    """
    mkdir -p output
    touch output/diamond_output.asn
    """
}
