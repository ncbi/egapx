#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params; to_map; shellSplit } from '../../utilities'
include { run_diamond_egap;} from '../../shared/diamond/main'


workflow diamond_worker {
    take:
        gnomon_prot_ids
        swiss_prot_ids
        gnomon_prot_asn
        swiss_prot_asn
        parameters  // Map : extra parameter and parameter update
    main:
        String diamond_blastp_params = merge_params('--sam-query-len --very-sensitive --unal 0 --comp-based-stats 0 --masking 0', parameters, 'diamond_blastp')
        String diamond_regular_params = merge_params('-ofmt seq-align-set -query-fmt seq-ids -subject-fmt seq-ids -output-prefix hits', parameters, 'diamond')
        String diamond_egap_params = '-blastp-args \''  + diamond_blastp_params + '\' '  + diamond_regular_params
        
        run_diamond_egap(gnomon_prot_ids, swiss_prot_ids, gnomon_prot_asn, swiss_prot_asn, diamond_egap_params)

    emit:
        alignments = run_diamond_egap.out
}
