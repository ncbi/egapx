#!/usr/bin/env nextflow
// main nextflow script for EGAPx execution
// route data to subworkflows

nextflow.enable.dsl=2

include { setup_genome; setup_proteins } from './setup/main'
include { sra_query } from './rnaseq_short/sra_qry/main'
include { fetch_sra_fasta } from './rnaseq_short/fetch_sra_fasta/main'
include { star_index } from './rnaseq_short/star_index/main'
include { star_wnode_simplified as star_simplified; star_wnode as star } from './rnaseq_short/star_wnode/main'
include { bam_strandedness } from './rnaseq_short/bam_strandedness/main'
include { bam_bin_and_sort } from './rnaseq_short/bam_bin_and_sort/main'
include { bam2asn } from './rnaseq_short/convert_from_bam/main'
include { rnaseq_collapse } from './rnaseq_short/rnaseq_collapse/main'
include { get_hmm_params } from './default/get_hmm_params/main'
include { chainer_wnode as chainer } from './gnomon/chainer_wnode/main'
include { gnomon_wnode } from './gnomon/gnomon_wnode/main'
include { prot_gnomon_prepare } from './gnomon/prot_gnomon_prepare/main'
include { annot_builder } from './default/annot_builder/main'
include { annotwriter } from './default/annotwriter/main'
include { miniprot } from './target_proteins/miniprot/main'
include { align_filter_sa } from './target_proteins/align_filter_sa/main'
include { best_aligned_prot } from './target_proteins/best_aligned_prot/main'
include { paf2asn } from './target_proteins/paf2asn/main'
include { run_align_sort} from './gnomon/align_sort_sa/main'

params.intermediate = false


workflow egapx {
    take:
        genome          // path to genome
        proteins        // path to proteins, optional
        // Alternative groups of parameters, one of them should be set
        // reads_query - SRA query in the form accepted by NCBI
        // reads_ids - list of SRA IDs
        // reads, reads_metadata - path to reads accompanied by metadata
        reads_query     // SRA query
        reads_ids       // list of SRA IDs
        reads           // path to reads
        reads_metadata  // path to reads metadata 13 tab-delimited fields, 1-st - SRA ID, 3-rd paired or unpaired, everything else - not used, but must be present
                        // 4, 5, 13 - numbers, 5 - non zero number
        organelles      // path to organelle list
        // Alternative parameters, one of them should be set
        // tax_id - NCBI tax id of the closest taxon to the genome
        // hmm_params - HMM parameters
        tax_id          // NCBI tax id of the closest taxon to the genome
        hmm_params      // HMM parameters
        //
        softmask        // softmask for GNOMON, optional
        task_params     // task parameters for every task
    main:
        def (scaffolds, gencoll_asn, unpacked_genome, genome_asn) = setup_genome(genome, organelles, task_params.get('setup', [:]))

        // Protein alignments
        def protein_alignments = []
        def unpacked_proteins
        def proteins_asn = []
        if (proteins) {
            // miniprot plane
            (unpacked_proteins, proteins_asn) = setup_proteins(proteins, task_params.get('setup', [:]))
            miniprot(unpacked_genome, unpacked_proteins, task_params.get('miniprot', [:]))
            paf2asn(genome_asn, proteins_asn, miniprot.out.miniprot_file, task_params.get('paf2asn', [:]))
            best_aligned_prot(genome_asn, proteins_asn, paf2asn.out.asn_file, gencoll_asn, task_params.get('best_aligned_prot', [:]))
            align_filter_sa(genome_asn, proteins_asn, best_aligned_prot.out.asn_file, task_params.get('align_filter_sa', [:]))
            protein_alignments = run_align_sort(genome_asn, proteins_asn,align_filter_sa.out.filtered_file, 
                "-k subject,subject_start,-subject_end,subject_strand,query,query_start,-query_end,query_strand,-num_ident,gap_count" )
        }

        // RNASeq short alignments
        def rnaseq_alignments = []
        // Satisfy quirks of Nextflow compiler
        def reads_query1 = reads_query
        def reads_ids1 = reads_ids
        // Conditional code on SRA reads source
        if (reads_query || reads_ids || reads) {
            def index = star_index(unpacked_genome, task_params.get('star_index', [:]))
            def ch_align, ch_align_index, sra_metadata, sra_run_list
            if (reads_query || reads_ids) {
                def query = reads_query1 ? reads_query1 : reads_ids1.join("[Accession] OR ") + "[Accession]"
                (sra_metadata, sra_run_list) = sra_query(query, task_params.get('sra_qry', [:]))
                def reads_fasta_pairs = fetch_sra_fasta(sra_run_list, task_params.get('fetch_sra_fasta', [:]))
                (ch_align, ch_align_index) = star(scaffolds, reads_fasta_pairs, genome_asn, index, task_params.get('star_wnode', [:]))
            } else {
                sra_metadata = reads_metadata
                (ch_align, ch_align_index) = star_simplified(scaffolds, reads, genome_asn, index, task_params.get('star_wnode', [:]))
            }
            //

            bam_strandedness(ch_align, ch_align_index, sra_metadata, task_params.get('bam_strandedness', [:]))
            def strandedness = bam_strandedness.out.strandedness
            
            // Run bam_bin_and_sort
            bam_bin_and_sort(ch_align, ch_align_index, unpacked_genome, organelles, task_params.get('bam_bin_and_sort', [:]))
            def bam_bins = bam_bin_and_sort.out.sorted

            // Run BAM2ASN
            bam2asn(bam_bins, strandedness, genome_asn, task_params.get('convert_from_bam', [:]))
            def asn_align = bam2asn.out.align.collect()
            def keylist = bam2asn.out.keylist.collect()

            rnaseq_collapse(genome_asn, keylist, asn_align, sra_metadata, 10, task_params.get('rnaseq_collapse', [:]))
            rnaseq_alignments = rnaseq_collapse.out.alignments
        }

        // Combine RNASeq and protein alignments

        def alignments
        if (proteins && (reads_query || reads_ids || reads)) [
            alignments = rnaseq_alignments.combine(protein_alignments)
        ] else if (proteins) {
            alignments = protein_alignments
        } else {
            alignments = rnaseq_alignments
        }

        // GNOMON

        def effective_hmm
        if (hmm_params) {
            effective_hmm = hmm_params
        } else {
            effective_hmm = get_hmm_params(tax_id, [:])
        }

        chainer(alignments, effective_hmm, /* evidence_denylist */ [], /* gap_fill_allowlist */ [], scaffolds, /* trusted_genes */ [], genome_asn, proteins_asn, task_params.get('chainer', [:]))

        gnomon_wnode(scaffolds, chainer.out.chains, chainer.out.chains_slices, effective_hmm, [], softmask, genome_asn, proteins_asn, task_params.get('gnomon', [:]))
        def models = gnomon_wnode.out.outputs

        // prot_gnomon_prepare(models, task_params.get('prot_gnomon_prepare', [:]))
        
        annot_builder(gencoll_asn, models, genome_asn, task_params.get('annot_builder', [:]))
        def accept_asn = annot_builder.out.accept_asn
        
        annotwriter(accept_asn, [:])
        annotwriter.out.annoted_file
    emit:
        out_files = annotwriter.out.annoted_file
        annot_builder_output = annot_builder.out.outputs
}
