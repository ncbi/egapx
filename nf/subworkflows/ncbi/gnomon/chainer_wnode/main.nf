#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { merge_params } from '../../utilities'
include { run_align_sort }  from '../../default/align_sort_sa/main.nf'


workflow chainer_wnode {
    take:
        sorted_alignments
        hmm_params
        evidence_denylist
        gap_fill_allowlist
        scaffolds
        trusted_genes    
        genome
        proteins
        parameters  // Map : extra parameter and parameter update
    main:
        String submit_chainer_params =  merge_params("-minimum-abut-margin 20 -separate-within-introns", parameters, 'submit_chainer')
        String chainer_wnode_params = merge_params("", parameters, 'chainer_wnode')
        String gpx_make_outputs_params =  merge_params("-default-output-name chains -slices-for affinity -sort-by affinity", parameters, 'gpx_make_outputs')
        
        def (jobs, lines_per_file) = generate_jobs(sorted_alignments, submit_chainer_params)
        def collected = run_chainer(jobs.flatten(), sorted_alignments, hmm_params, evidence_denylist, gap_fill_allowlist, scaffolds, trusted_genes, genome, proteins, lines_per_file, chainer_wnode_params) | collect

        run_gpx_make_outputs(collected, gpx_make_outputs_params)
    emit:
        chains = run_gpx_make_outputs.out.chains
        chains_slices = run_gpx_make_outputs.out.chains_slices
        evidence = run_gpx_make_outputs.out.evidence
        evidence_slices = run_gpx_make_outputs.out.evidence_slices
        all = run_gpx_make_outputs.out.all
}


process generate_jobs {
    label 'gpx_submitter'
    label 'small_mem'
    input:
        path sort_aligns
        val params
    output:
        path "job.*"
        env lines_per_file
    script:
        njobs=task.ext.split_jobs
    """
    #!/usr/bin/env bash
    # generate_jobs $sort_aligns $params -output chains -output-slices chains_slices -output-evidence evidence -output-evidence-slices evidence_slices
    submit_chainer $params -asn $sort_aligns -o jobs
    total_lines=\$(wc -l <jobs)
    (( lines_per_file = (total_lines + ${njobs} - 1) / ${njobs} ))
    echo total_lines=\$total_lines, lines_per_file=\$lines_per_file
    ####split -l\$lines_per_file jobs job. -da 3
    # Use round robin to distribute jobs across nodes more evenly
    if [ \$total_lines -lt $njobs ]; then
        effective_njobs=\$total_lines
    else
        effective_njobs=$njobs
    fi
    split -nr/\$effective_njobs jobs job. -da 3
    """
    stub:
    njobs=task.ext.split_jobs
    """
    for i in {1..$njobs}; do
        echo "<job query =\\\"lcl|${sort_aligns}:\${i}-\${i}\\\"></job>" >> jobs
    done
    split -nr/$njobs jobs job. -da 3
    lines_per_file=10
    """
}


process run_chainer {   
    label 'long_job'
    label 'multi_node'
    label 'med_mem'
    input:
        path job
        path alignments
        path hmm_params
        path evidence_denylist
        path gap_fill_allowlist
        path scaffolds
        path trusted_genes
        path genome, stageAs: 'indexed/*'
        path proteins_asn, stageAs: 'indexed/*'
        val lines_per_file
        val params
    output:
        path "output/*"
    script:
    """
    echo "${evidence_denylist.join('\n')}" > evidence_denylist.mft
    echo "${gap_fill_allowlist.join('\n')}" > gap_fill_allowlist.mft
    echo "${scaffolds.join('\n')}" > scaffolds.mft
    echo "${trusted_genes.join('\n')}" > trusted_genes.mft
    # HACK: derive start_job_id from job file extension
    filename=\$(basename -- "$job")
    extension="\${filename##*.}"
    (( start_job_id = ((10#\$extension) * $lines_per_file) + 1 ))
    
    mkdir -p tmp
    # make the local LDS of the genomic and protein (if present) sequences
    ##lds2_indexer -source indexed -db tmp/LDS2
    mkdir -p tmp/asncache/
    auto_prime_cache.py -cache tmp/asncache/ -i ${genome} -oseq-ids spidsg -split-sequences
    auto_prime_cache.py -cache tmp/asncache/ -i ${proteins_asn} -oseq-ids spidsp -split-sequences


    mkdir -p tmp/interim
    chainer_wnode $params -start-job-id \$start_job_id  -workers ${task.ext.threads} -input-jobs ${job} -O tmp/interim -nogenbank -asn-cache tmp/asncache/ -evidence-denylist-manifest evidence_denylist.mft -gap-fill-allowlist-manifest gap_fill_allowlist.mft -param ${hmm_params} -scaffolds-manifest scaffolds.mft -trusted-genes-manifest trusted_genes.mft
    mkdir -p output
    cat tmp/interim/* > output/chainer_wnode.${task.index}.gpx-job.asnb
    rm -rf tmp
    """

    stub:
    """
    mkdir -p output
    touch output/chainer_wnode.${task.index}.gpx-job.asnb
    """
}


process run_gpx_make_outputs {
    label 'single_cpu'
    label 'small_mem'
    input:
        path files, stageAs: "gpx_inputs/*"
        val params
    output:
        path "output/chains.*.out.gz", emit: 'chains'
        path "output/chains.*.out.gz.slices", emit: 'chains_slices'
        path "output/evidence.*.out.gz", emit: 'evidence', optional: true
        path "output/evidence.*.out.gz.slices", emit: 'evidence_slices', optional: true
        path "output/*", emit: 'all'
    script:
    """
    ls -1 gpx_inputs/* > gpx_inputs.mft
    mkdir -p output
    gpx_make_outputs $params -input-manifest gpx_inputs.mft -output output/@.#.out.gz -output-manifest output/@.mft -slices-manifest output/@_slices.mft -num-partitions ${task.ext.split_jobs}
    """
    stub:
    """
    mkdir -p output
    echo ${files}
    for i in {1..${task.ext.split_jobs} }; do
        touch output/chains.\$i.out.gz
        touch output/chains.\$i.out.gz.slices
        touch output/evidence.\$i.out.gz
        touch output/evidence.\$i.out.gz.slices
    done
    """
}
