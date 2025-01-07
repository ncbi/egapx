#!/usr/bin/env nextflow
// main nextflow script for EGAPx ui execution
// prepare data channels and call main subworkflow

nextflow.enable.dsl=2

include { egapx } from './subworkflows/ncbi/main'

params.verbose = false


process export {
    publishDir "${params.output}", mode: 'copy', saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    input:
        path genomic_gff
        path genomic_gtf
        path genomic_fasta
        path transcript_fasta
        path cds_fasta
        path proteins_fasta
        path annot_builder_output, stageAs: 'annot_builder_output/*'
        path validated, stageAs: 'validated/*'
        path stats, stageAs: 'stats/*'
        path annotated_genome_asn
        path annotation_data_comment
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
    def proteins_trusted = input_params.get('proteins_trusted', [])
    def reads_query = input_params.get('reads_query', [])
    def reads_ids = input_params.get('reads_ids', [])
    def reads = input_params.get('reads', [])
    def reads_metadata = input_params.get('reads_metadata', [])
    def organelles = input_params.get('organelles', []) ?: []
    def tax_id = input_params.get('taxid', [])
    def symbol_format_class = input_params.get('symbol_format_class', [])
    def hmm_params = input_params.get('hmm', []) ?: []
    def train_hmm = input_params.get('train_hmm', [])
    def softmask = input_params.get('softmask', []) ?: []
    def max_intron = input_params.get('max_intron', [])
    def genome_size_threshold = input_params.get('genome_size_threshold', [])
    def ortho_files = input_params.get('ortho', []) ?: []
    def rnaseq_alignments = input_params.get('rnaseq_alignments', []) ?: []
    def protein_alignments = input_params.get('protein_alignments', []) ?: []
    def reference_sets = input_params.get('reference_sets', []) ?: []
    def prot_denylist = input_params.get('prot_denylist', []) ?: []
    def task_params = params.get('tasks', [:])
    if (params.verbose) {
        println("input params:\ngenome ${genome}")
        println("proteins ${proteins}")
        println("proteins_trusted ${proteins_trusted}")
        println("reads_query ${reads_query}")
        println("reads_ids ${reads_ids}")
        println("reads ${reads}")
        println("reads_metadata ${reads_metadata}")
        println("organelles ${organelles}")
        println("tax_id ${tax_id}")
        println("symbol_format_class ${symbol_format_class}")
        println("hmm_params ${hmm_params}")
        println("train_hmm ${train_hmm}")
        println("softmask ${softmask}")
        println("max_intron ${max_intron}")
        println("genome_size_threshold ${genome_size_threshold}")
        println("ortho_files ${ortho_files}")
        println("rnaseq_alignments ${rnaseq_alignments}")
        println("protein_alignments ${protein_alignments}")
        println("reference_sets ${reference_sets}")
        println("prot_denylist ${prot_denylist}")
        // Keep it last as it is large
        println("task_params ${task_params}")
    }
    
    egapx(genome, proteins, proteins_trusted, reads_query, reads_ids, reads, reads_metadata, organelles, tax_id, symbol_format_class, hmm_params, train_hmm, softmask, max_intron, genome_size_threshold, ortho_files, reference_sets, prot_denylist, task_params)
    export(egapx.out.out_ggff, 
           egapx.out.out_ggtf,
           egapx.out.out_gfa,
           egapx.out.out_rna_fa,
           egapx.out.out_cds_fa,
           egapx.out.out_prot_fa,
           egapx.out.annot_builder_output,
           egapx.out.validated,
           egapx.out.stats,
           egapx.out.annotated_genome_asn,
           egapx.out.annotation_data_comment)
}
