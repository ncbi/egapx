#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow gnomon_wnode {
    take:
        scaffolds
        chains
        chains_slices
        hmm_params
        softmask_lds2
        softmask_lds2_source
        genome
        proteins
        parameters  // Map : extra parameter and parameter update
    main:
        String gpx_qsubmit_params =  merge_params("", parameters, 'gpx_qsubmit')
        String annot_params =  merge_params("-margin 1000 -mincont 1000 -minlen 225 -mpp 10.0 -ncsp 25 -window 200000 -nonconsens -open", parameters, 'annot_wnode')
        String gpx_qdump_params =  merge_params("-slices-for affinity -sort-by affinity", parameters, 'gpx_qdump')

        def (jobs, lines_per_file) = gpx_qsubmit(scaffolds, chains, chains_slices, gpx_qsubmit_params)
        def annot_files = annot(jobs.flatten(), chains, hmm_params, softmask_lds2, softmask_lds2_source, genome, proteins, lines_per_file, annot_params)
        gpx_qdump(annot_files.collect(), gpx_qdump_params)
    emit:
        outputs = gpx_qdump.out.outputs
}


process gpx_qsubmit {
    input:
        path scaffolds
        path chains
        path chains_slices
        val params
    output:
        path "job.*"
        env lines_per_file
    script:
        njobs=16
    """
    echo $scaffolds | tr ' ' '\\n' > scaffolds.mft
    for file in $chains_slices; do
        echo \$file >> chains_slices.mft
        # remove path from the first line of this file
        sed -i -e '1s/\\(.*\\)\\/\\(.*\\)\$/\\2/' \$file
    done
    gpx_qsubmit $params -ids-manifest scaffolds.mft -slices-manifest chains_slices.mft -o jobs
    total_lines=\$(wc -l <jobs)
    (( lines_per_file = (total_lines + ${njobs} - 1) / ${njobs} ))
    echo total_lines=\$total_lines, lines_per_file=\$lines_per_file
    # split -l\$lines_per_file jobs job. -da 3
    # Use round robin to distribute jobs across nodes more evenly
    if [ \$total_lines -lt $njobs ]; then
        effective_njobs=\$total_lines
    else
        effective_njobs=$njobs
    fi
    split -nr/\$effective_njobs jobs job. -da 3
    """
}


process annot {
    input:
        path jobs
        path chains // used for staging chain files, referred from jobs
        path hmm_params
        path softmask_lds2
        path softmask
        path genome, stageAs: 'indexed/*'
        path proteins_asn, stageAs: 'indexed/*'
        val lines_per_file
        val params
    output:
        path "outdir/*"
    script:
    """
    njobs=`wc -l <$jobs`
    if [ \$njobs -lt 16 ]; then
        threads=\$njobs
    else
        threads=16
    fi

    lds2=indexed_lds
    if [ -n "$softmask_lds2" ]; then
        # patch LDS2 to point to the source
        files=\$(sqlite3 $softmask_lds2 -cmd "SELECT file_name FROM file" ".exit")
        for f in \$files; do
            base=\$(basename \$f)
            sqlite3 $softmask_lds2 -cmd "UPDATE file SET file_name = '\$base' WHERE file_name = '\$f'" ".exit"
        done
        lds2+=",$softmask_lds2"
    elif [ -n "$softmask" ]; then
        mkdir sm_src
        mv $softmask ./sm_src/
        lds2_indexer -source ./sm_src/ -db softmask_lds2
        lds2+=",softmask_lds2"
    fi    

    filename=\$(basename -- "$jobs")
    extension="\${filename##*.}"
    (( start_job_id = ((10#\$extension) * $lines_per_file) + 1 ))

    # make the local LDS of the genomic fasta
    lds2_indexer -source indexed -db indexed_lds

    # When running multiple jobs on the cluster there is a chance that
    # several jobs will run on the same node and thus generate files
    # with the same filename. We need to avoid that to be able to stage
    # the output files for gpx_make_outputs. We add the job file numeric
    # extension as a prefix to the filename.
    mkdir interim
    annot_wnode $params -nogenbank -lds2 \$lds2  -start-job-id \$start_job_id -workers \$threads -input-jobs $jobs -param $hmm_params -O interim || true
    mkdir outdir
    for f in interim/*; do
        if [ -f \$f ]; then
            mv \$f outdir/\${extension}_\$(basename \$f)
        fi
    done
    """
}


process gpx_qdump {
    input:
        path files, stageAs: "inputs/*"
        val params
    output:
        path "*.out", emit: "outputs"
    script:
    """
    gpx_qdump $params -input-path inputs -output gnomon_wnode.out
    """
}
