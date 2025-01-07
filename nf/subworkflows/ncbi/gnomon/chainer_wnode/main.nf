#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include { merge_params } from '../../utilities'
include { run_align_sort }  from '../../default/align_sort_sa/main.nf'

split_count=16


workflow chainer_wnode {
    take:
        alignments
        hmm_params
        evidence_denylist
        gap_fill_allowlist
        scaffolds
        trusted_genes    
        genome
        proteins
        parameters  // Map : extra parameter and parameter update
    main:
        String input_sorting = parameters.get('input_aligns_sort', '')
        def sort_aligns = alignments
        if (!input_sorting.contains("presorted")) {
            String align_sort_params = ""
            if (input_sorting.contains("merge_only")) {
                align_sort_params = "-merge"
            }
            align_sort_params += " -ifmt seq-align -compression none -k subject,subject_start,-subject_end "
            // print(align_sort_params)
            sort_aligns = run_align_sort([], [], alignments, align_sort_params).collect()
            //sort_aligns = align_sort(alignments, align_sort_params)
        }
        String submit_chainer_params =  merge_params("-minimum-abut-margin 20 -separate-within-introns", parameters, 'submit_chainer')
        String chainer_wnode_params = merge_params("", parameters, 'chainer_wnode')
        String gpx_make_outputs_params =  merge_params("-default-output-name chains -slices-for affinity -sort-by affinity", parameters, 'gpx_make_outputs')
        
        def (jobs, lines_per_file) = generate_jobs(sort_aligns, submit_chainer_params)
        def collected = run_chainer(jobs.flatten(), sort_aligns, hmm_params, evidence_denylist, gap_fill_allowlist, scaffolds, trusted_genes, genome, proteins, lines_per_file, chainer_wnode_params) | collect

        run_gpx_make_outputs(collected, gpx_make_outputs_params)
    emit:
        chains = run_gpx_make_outputs.out.chains
        chains_slices = run_gpx_make_outputs.out.chains_slices
        evidence = run_gpx_make_outputs.out.evidence
        evidence_slices = run_gpx_make_outputs.out.evidence_slices
        all = run_gpx_make_outputs.out.all
}


process generate_jobs {
    input:
        path sort_aligns
        val params
    output:
        path "job.*"
        env lines_per_file
    script:
        njobs=split_count
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
    """
    for i in {1..$split_count}; do
        echo "<job query =\\\"lcl|${sort_aligns}:\${i}-\${i}\\\"></job>" >> jobs
    done
    split -nr/$split_count jobs job. -da 3
    lines_per_file=10
    """
}


process run_chainer {
    label 'long_job'
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
        job_num = job.toString().tokenize('.').last().toInteger()
    """
    echo "${evidence_denylist.join('\n')}" > evidence_denylist.mft
    echo "${gap_fill_allowlist.join('\n')}" > gap_fill_allowlist.mft
    echo "${scaffolds.join('\n')}" > scaffolds.mft
    echo "${trusted_genes.join('\n')}" > trusted_genes.mft
    # HACK: derive start_job_id from job file extension
    filename=\$(basename -- "$job")
    extension="\${filename##*.}"
    (( start_job_id = ((10#\$extension) * $lines_per_file) + 1 ))
    
    # make the local LDS of the genomic and protein (if present) sequences
    lds2_indexer -source indexed -db LDS2

    # When running multiple jobs on the cluster there is a chance that
    # several jobs will run on the same node and thus generate files
    # with the same filename. We need to avoid that to be able to stage
    # the output files for gpx_make_outputs. We add the job file numeric
    # extension as a prefix to the filename.
    mkdir -p interim
    chainer_wnode $params -start-job-id \$start_job_id  -workers 32 -input-jobs ${job} -O interim -nogenbank -lds2 LDS2 -evidence-denylist-manifest evidence_denylist.mft -gap-fill-allowlist-manifest gap_fill_allowlist.mft -param ${hmm_params} -scaffolds-manifest scaffolds.mft -trusted-genes-manifest trusted_genes.mft
    mkdir -p output
    for f in interim/*; do
        if [ -f \$f ]; then
            mv \$f output/\${extension}_\$(basename \$f)
        fi
    done
    """

    stub:
        job_num = job.toString().tokenize('.').last().toInteger()
    """
    mkdir -p output
    touch output/sample_chainer_wnode.${job_num}.out
    """
}


process run_gpx_make_outputs {
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
    gpx_make_outputs $params -input-manifest gpx_inputs.mft -output output/@.#.out.gz -output-manifest output/@.mft -slices-manifest output/@_slices.mft -num-partitions $split_count
    """
    stub:
    """
    mkdir -p output
    echo ${files}
    for i in {1..$split_count}; do
        touch output/chains.\$i.out.gz
        touch output/chains.\$i.out.gz.slices
        touch output/evidence.\$i.out.gz
        touch output/evidence.\$i.out.gz.slices
    done
    """
}
