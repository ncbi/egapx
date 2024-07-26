#!/usr/bin/env nextflow
// gnomon plane workflow
// route data to tasks

nextflow.enable.dsl=2

include { chainer_wnode as chainer } from './chainer_wnode/main'
include { gnomon_wnode } from './gnomon_wnode/main'
include { prot_gnomon_prepare } from './prot_gnomon_prepare/main'
include { gnomon_training_iterations } from '../gnomon-training-iteration/main'

include { diamond_worker} from './diamond/main'
include { best_protein_hits } from './protein_filter/main'
include { gnomon_biotype} from './gnomon_biotype/main'
include { fetch_swiss_prot_asn; get_swiss_prot_ids } from '../shared/diamond/main'
include { diamond_orthology } from '../orthology/diamond_orthology/main'
include { locus_link } from './locus_link/main'


params.intermediate = false

workflow gnomon_plane {
    take:
        genome_asn
        scaffolds
        gencoll_asn
        proteins_asn
        alignments // list of all relevent input alignments

        // Alternative parameters, one of them should be set
        // tax_id - NCBI tax id of the closest taxon to the genome
        // hmm_params - HMM parameters
        tax_id          // NCBI tax id of the closest taxon to the genome
        hmm_params      // HMM parameters
        hmm_taxid       // NCBI tax id of the taxon of the HMM
        //
        softmask        // softmask for GNOMON, optional
        max_intron      // max intron length
        task_params     // task parameters for every task
    main:
        // GNOMON
        def effective_hmm
        if (tax_id == hmm_taxid) {
            effective_hmm = hmm_params
        } else {
            effective_hmm = gnomon_training_iterations(hmm_params, genome_asn, proteins_asn, alignments, /* evidence_denylist */ [], /* gap_fill_allowlist */ [],
                /* trusted_genes */ [], scaffolds, softmask,
                softmask, scaffolds,
                max_intron,
                task_params)
        }

        chainer(alignments, effective_hmm, /* evidence_denylist */ [], /* gap_fill_allowlist */ [], scaffolds, /* trusted_genes */ [], genome_asn, proteins_asn, task_params.get('chainer', [:]))

        def gn_models = []
        gnomon_wnode(scaffolds, chainer.out.chains, chainer.out.chains_slices, effective_hmm, [], softmask, genome_asn, proteins_asn, task_params.get('gnomon', [:]))
       
    emit:
        gnomon_models = gnomon_wnode.out.outputs
        // trained_hmm = effective_hmm
}



workflow post_gnomon_plane {
    take:
        gnomon_models
        gencoll_asn
        orthologs


        // Alternative parameters, one of them should be set
        // tax_id - NCBI tax id of the closest taxon to the genome
        // hmm_params - HMM parameters
        tax_id          // NCBI tax id of the closest taxon to the genome
        task_params     // task parameters for every task
    main:
        // Post GNOMON
        // might come its own plane
        def swiss_prot_asn = fetch_swiss_prot_asn()
        def swiss_prot_ids = get_swiss_prot_ids(swiss_prot_asn)

        prot_gnomon_prepare(gnomon_models, task_params.get('prot_gnomon_prepare', [:]))
        // Seed Protein-Model Hits
        diamond_worker(prot_gnomon_prepare.out.prot_ids, swiss_prot_ids, gnomon_models, swiss_prot_asn, task_params.get('diamond', [:]))
        best_protein_hits(gnomon_models, swiss_prot_asn,  diamond_worker.out.alignments , task_params.get('protein_filter', [:]))

        gnomon_biotype([] /*models*/,/*splices_file  -- constant*/ [],  /*denylist -- constant*/ [], gencoll_asn, swiss_prot_asn, gnomon_models, diamond_worker.out.alignments,task_params.get('gnomon_biotype', [:]))
        locus_link(/*best_refseq_prot_hit  -- best protein hits from refseq plane*/ [], orthologs, []  /*annot_builder.out.annot_files*/, 
                gencoll_asn, gnomon_models, best_protein_hits.out.alignments , /*track_loci*/ [], /*comparisons*/ [],  /*curr_prev_compare*/ [], 
                gnomon_biotype.out.biotypes, /*lxr_data*/ [], swiss_prot_asn, /*name_from_ortholog */ [],  task_params.get('locus_link', [:]))

       
    emit:
        locus = locus_link.out.locus
}
