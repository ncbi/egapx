#!/usr/bin/env nextflow
// main nextflow script for EGAPx ui execution
// prepare data channels and call main subworkflow

nextflow.enable.dsl=2

include { egapx } from './subworkflows/ncbi/main'
include { only_gnomon } from './subworkflows/ncbi/only_gnomon'

params.verbose = false


process export {
    publishDir "${params.output}", mode: 'copy', saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    input:
        path out_files
        path annot_builder_output, stageAs: 'annot_builder_output/*'
        // path locus
    output:
        path "*", includeInputs: true
    script:
    """
    echo "export script"
    """
    stub:
    """
    echo "export stub"
    """
}


workflow {
    // Parse input parameters
    def input_params = params.get('input', [:])
    def genome = input_params.get('genome', [])
    def proteins = input_params.get('proteins', [])
    def reads_query = input_params.get('reads_query', [])
    def reads_ids = input_params.get('reads_ids', [])
    def reads = input_params.get('reads', [])
    def reads_metadata = input_params.get('reads_metadata', [])
    def organelles = input_params.get('organelles', []) ?: []
    def tax_id = input_params.get('taxid', [])
    def hmm_params = input_params.get('hmm', []) ?: []
    def hmm_taxid = input_params.get('hmm_taxid', []) ?: []
    def softmask = input_params.get('softmask', []) ?: []
    def max_intron = input_params.get('max_intron', [])
    def genome_size_threshold = input_params.get('genome_size_threshold', [])
    def rnaseq_alignments = input_params.get('rnaseq_alignments', []) ?: []
    def protein_alignments = input_params.get('protein_alignments', []) ?: []
    def task_params = params.get('tasks', [:])
    def func_name = params.get('func_name', '')
    if (params.verbose) {
        println("input params:\ngenome ${genome}")
        println("proteins ${proteins}")
        println("reads_query ${reads_query}")
        println("reads_ids ${reads_ids}")
        println("reads ${reads}")
        println("reads_metadata ${reads_metadata}")
        println("organelles ${organelles}")
        println("tax_id ${tax_id}")
        println("hmm_params ${hmm_params}")
        println("hmm_taxid ${hmm_taxid}")
        println("softmask ${softmask}")
        println("max_intron ${max_intron}")
        println("genome_size_threshold ${genome_size_threshold}")
        println("rnaseq_alignments ${rnaseq_alignments}")
        println("protein_alignments ${protein_alignments}")
        println("func_name ${func_name}")
        // Keep it last as it is large
        println("task_params ${task_params}")
    }
    
    if(func_name == 'only_gnomon') {
        if (params.verbose) {
            print('in gnomon block')
        }
        only_gnomon(genome, proteins, rnaseq_alignments, protein_alignments, organelles, tax_id, hmm_params, hmm_taxid, softmask, task_params)
        export(only_gnomon.out.out_files, only_gnomon.out.evidence)
    }
    else {
        if (params.verbose) {
            print('in egapx block')
        }
        egapx(genome, proteins, reads_query, reads_ids, reads, reads_metadata, organelles, tax_id, hmm_params, hmm_taxid, softmask, max_intron, genome_size_threshold, task_params)
        // export(egapx.out.out_files, egapx.out.annot_builder_output, egapx.out.locus)
        export(egapx.out.out_files, egapx.out.annot_builder_output)
    }
}
