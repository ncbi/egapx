#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow rnaseq_collapse {
    take:
        genome            // path: file of genomic sequences, FASTA or ASN
        scaffold_list     // path: file with list of scaffolds
        alignments        // value channel path: align files
        sra_metadata_list // path: file with sra metadata
        max_jobs          // val: max number of jobs
        parameters        // Map: extra parameter and parameter update
    main:
        String create_jobs_params =  merge_params("-alignments-per-job 50000 -min-range 100000", parameters, 'rnaseq_collapse_create_jobs')
        String rnaseq_collapse_params = merge_params("-backlog 1 -max-jobs 1 -support-non-sra", parameters, 'rnaseq_collapse')
        String gpx_make_outputs_params =  merge_params("-default-output-name align -slices-for affinity -sort-by job-id -unzip align", parameters, 'gpx_make_outputs')

        def (jobs, lines_per_file) = generate_jobs(genome, scaffold_list, max_jobs, create_jobs_params)
        def collected = run_rnaseq_collapse(genome, scaffold_list, alignments, sra_metadata_list, jobs.flatten(), lines_per_file, rnaseq_collapse_params) | collect
        run_gpx_make_outputs(collected, gpx_make_outputs_params)
    emit:
        alignments = run_gpx_make_outputs.out.alignments
        alignment_slices = run_gpx_make_outputs.out.alignment_slices
}


process generate_jobs {
    input:
        path genome, stageAs: 'genome/*'
        path scaffold_list
        val  njobs
        val  params
    output:
        path "job.*"
        env lines_per_file
    script:
    """
    echo "${scaffold_list.join('\n')}" > scaffold_list.mft
    
    # make the local LDS of the genomic sequences
    lds2_indexer -source ./genome -db ./genome_lds  

    rnaseq_collapse_create_jobs $params -nogenbank -lds2 ./genome_lds -scaffold-list scaffold_list.mft > jobs
    total_lines=\$(wc -l <jobs)
    (( lines_per_file = (total_lines + ${njobs} - 1) / ${njobs} ))
    split -l\$lines_per_file jobs job. -da 3
    """

    stub:
    """
    for i in {1..$njobs}; do
        echo "<job query =\\\"lcl|SOME_GENOME_ID:\${i}-\${i}\\\"></job>" >> jobs
    done
    split -l1 jobs job. -da 3
    lines_per_file=1
    """
}


process run_rnaseq_collapse {
    input:
        path genome, stageAs: 'genome/*'
        path scaffold_list
        path in_align
        path sra_metadata_list
        path job
        val  lines_per_file
        val  params
    output:
        path "output/*"
    script:
    """
    njobs=`wc -l <$job`
    if [ \$njobs -lt 16 ]; then
        threads=\$njobs
    else
        threads=16
    fi

    echo "${scaffold_list.join('\n')}" > scaffold_list.mft
    echo "${in_align.join('\n')}" > align.mft
    echo "${sra_metadata_list.join('\n')}" > metadata.mft

    # HACK: derive start_job_id from job file extension
    filename=\$(basename -- "$job")
    extension="\${filename##*.}"
    # NB: for successful gather phase all job id should be unique,
    # so we must supply start_job_id.
    (( start_job_id = ((10#\$extension) * $lines_per_file) + 1 ))
  
    # make the local LDS of the genomic sequences
    lds2_indexer -source ./genome -db ./genome_lds  
  
    # When running multiple jobs on the cluster there is a chance that
    # several jobs will run on the same node and thus generate files
    # with the same filename. We need to avoid that to be able to stage
    # the output files for gpx_make_outputs. We add the job file numeric
    # extension as a prefix to the filename.
    mkdir -p interim
    rnaseq_collapse $params -O interim -nogenbank -lds2 ./genome_lds -sorted-vols align.mft -scaffold-list scaffold_list.mft -sra-metadata-manifest metadata.mft -start-job-id \$start_job_id -input-jobs $job -workers \$threads
    mkdir -p output
    for f in interim/*; do
        if [ -f \$f ]; then
            mv \$f output/\${extension}_\$(basename \$f)
        fi
    done
    """
    
    stub:
    """
    filename=\$(basename -- "$job")
    extension="\${filename##*.}"
    mkdir -p output
    touch output/rnaseq_collapse_wnode.\${extension}.out
    """

}


process run_gpx_make_outputs {
    input:
        path files, stageAs: 'input/*'
        val  params
    output:
        path "output/align.*.out", emit: 'alignments'
        path "output/align.*.out.slices", emit: 'alignment_slices'
    script:
    """
    mkdir -p output
    gpx_make_outputs $params -input-path input -output output/@.#.out -output-manifest output/@.mft -slices-manifest output/@_slices.mft -num-partitions 1
    """

    stub:
    """
    mkdir -p output
    echo ${files}
    for i in {1..10}; do
        touch output/align.\$i.out
        touch output/align.\$i.out.slices
    done

    """
}
