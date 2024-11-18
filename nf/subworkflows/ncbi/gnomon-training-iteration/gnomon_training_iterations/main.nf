#!/usr/bin/env nextflow
nextflow.enable.dsl=2
//nextflow.preview.recursion=true


include { gnomon_training_iteration; gnomon_training_iteration as gnomon_training_iteration2; 
            gnomon_training_iteration as gnomon_training_iteration3; gnomon_training_iteration as gnomon_training_iteration4 } from './../utilities'
            

workflow gnomon_training_iterations {
    take:
        models_file
        genome_asn
        proteins_asn
        chainer_alignments
        chainer_evidence_denylist
        chainer_gap_fill_allowlist
        chainer_trusted_genes
        chainer_scaffolds
        gnomon_softmask
        gnomon_scaffolds
        max_intron
        parameters
    main:
    gnomon_training_iteration(models_file, genome_asn, proteins_asn ,chainer_alignments,chainer_evidence_denylist,chainer_gap_fill_allowlist,
               chainer_trusted_genes, chainer_scaffolds, 
               gnomon_softmask, gnomon_scaffolds, max_intron, parameters)
    gnomon_training_iteration2(gnomon_training_iteration.out.hmm_params_file, genome_asn, proteins_asn ,chainer_alignments,
               chainer_evidence_denylist,chainer_gap_fill_allowlist, chainer_trusted_genes, chainer_scaffolds, 
               gnomon_softmask, gnomon_scaffolds, max_intron, parameters)
    gnomon_training_iteration3(gnomon_training_iteration2.out.hmm_params_file, genome_asn, proteins_asn ,chainer_alignments,
               chainer_evidence_denylist,chainer_gap_fill_allowlist, chainer_trusted_genes, chainer_scaffolds, 
               gnomon_softmask, gnomon_scaffolds, max_intron, parameters)
    gnomon_training_iteration4(gnomon_training_iteration3.out.hmm_params_file, genome_asn, proteins_asn ,chainer_alignments,
               chainer_evidence_denylist,chainer_gap_fill_allowlist, chainer_trusted_genes, chainer_scaffolds, 
               gnomon_softmask, gnomon_scaffolds, max_intron, parameters)

    emit:
        hmm_params_file = gnomon_training_iteration4.out.hmm_params_file
}

workflow gnomon_no_training {
    take:
        hmm
     main:
         process_no_training(hmm)
    emit:
        hmm_params_file = process_no_training.out.file
}

process process_no_training
{
  input:
  val   hmm
  output:
    path "*.params" , emit: 'file' 
  script:
  """
   curl -O  ${hmm} 
  """
  stub:
  """
  touch hmm.params
  """
}


/*

//experimental (preview.recursion)
// to be used in the future
workflow gnomon_training_iterations {
    take:
        models_file
        genome_asn
        proteins_asn
        chainer_alignments
        chainer_evidence_denylist
        chainer_gap_fill_allowlist
        chainer_trusted_genes
        chainer_scaffolds
        gnomon_softmask_lds2
        gnomon_softmask_lds2_source
        gnomon_scaffolds
        max_intron
        parameters
    main:
        gnomon_training_iteration
            .recurse(models_file, genome_asn, proteins_asn ,chainer_alignments,chainer_evidence_denylist,chainer_gap_fill_allowlist, 
               chainer_trusted_genes,  chainer_scaffolds, gnomon_softmask_lds2,
               gnomon_softmask_lds2_source, gnomon_scaffolds, max_intron, parameters)
            .times(4)
    emit:
        hmm_params_file = gnomon_training_iteration.out.hmm_params_file
}

*/
