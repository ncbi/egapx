#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow minimap2 {
    take:
        genome_fasta
        genome_index // minimap2 index
        gencoll
        reads       // channel: FASTA file pairs generated from SRA reads, see e.g., https://www.nextflow.io/docs/latest/channel.html#fromfilepairs
        max_intron  // max intron length
        parameters  // Map : extra parameter and parameter update
    main:
        String gpx_qsubmit_params = merge_params("", parameters, "gpx_qsubmit")
        def (jobs, read_files, ids, reads_LDS, jobs_per_file) = gpx_qsubmit(genome_index, reads, gpx_qsubmit_params)
        def aux_data = read_files.join(ids).join(reads_LDS).join(jobs_per_file)
        // aux_data.view()
        expanded_jobs = jobs.transpose().combine(aux_data, by: 0)
        // expanded_jobs.view()
        running_sum = expanded_jobs.reduce( [1] ) { acc, v ->
            // println v
            acc.plus([acc[-1] + v[5].toInteger()])
        }
        new_jobs = expanded_jobs.merge(running_sum.flatten()) { a, b -> [a[0], a[1], a[2], a[3], a[4], b] }
        // new_jobs.view()
        alignments = minimap2_wnode(genome_fasta, genome_index, gencoll, new_jobs, max_intron, parameters)
        gpx_qdump(alignments.collect())
    emit:
        alignments = gpx_qdump.out.alignments
}


process gpx_qsubmit {
    input:
        path genome_index
        tuple val(sampleID), path(fasta_rna_file, stageAs: "reads/*")
        val parameters
    output:
        tuple val(sampleID), path("jobs/*")
        tuple val(sampleID), path(fasta_rna_file)
        tuple val(sampleID), path("${sampleID}.ids")
        tuple val(sampleID), path("${sampleID}.LDS")
        tuple val(sampleID), stdout
    script:
    """
    #!/usr/bin/env python3
    import os
    import math
    import shlex
    import subprocess
    import sqlite3

    def read_parameter(parts, name, dflt):
        p = dflt
        try:
            i = parts.index(name)
            if i + 1 < len(parts):
                try:
                    p = int(parts[i+1])
                except ValueError:
                    pass
        except ValueError:
            pass
        return p

    os.makedirs("jobs", exist_ok=True)
    idfile = "${sampleID}.ids"
    LDSfile = "${sampleID}.LDS"
    subject = "${genome_index}"
    # Read batch_size parameter
    parts = shlex.split("${parameters}")
    batch_size = read_parameter(parts, "-batch-size", 1000)
    # As we call qsubmit per sample run, -max-execs limits the number of
    # generated job files per sample, not total
    max_execs = read_parameter(parts, "-max-execs", 100)
    # Make LDS2 index of reads
    subprocess.run(["lds2_indexer", "-db", LDSfile, "-source", "reads/"])
    # Make ids from LDS2 index, fix file_name to relative to 'reads'
    with sqlite3.connect(f"file:{LDSfile}?mode=rw", uri=True) as conn:
        cur = conn.cursor()
        # ids
        res = cur.execute("SELECT txt_id FROM seq_id WHERE orig=1 AND int_id IS NULL")
        ids = res.fetchall()
        with open(idfile, "wt") as ids_file:
            ids_file.write(os.linesep.join(map(lambda x: x[0], ids)))
        # fix file_name
        res = cur.execute("SELECT file_name FROM file")
        files = res.fetchall()
        for f in files:
            new_name = "reads/" + os.path.basename(f[0])
            cur.execute("UPDATE file SET file_name=? WHERE file_name=?", (new_name, f[0]))
        conn.commit()
    num_reads = len(ids)
    batches_per_jobfile = math.ceil(num_reads / batch_size / max_execs)
    # redundant read, but we use tested code here and provide for possible
    # comments in ids file
    with open(idfile, "rb") as f:
        cntnt = f.read().decode('utf-8')
    job_num = 0
    offs = 0
    chunk_len = 0
    n_ids = 0
    jobs_for_file = []
    for line in cntnt.split(os.linesep):
        if not line or line.startswith("#"):
            offs += len(line) + 1
            continue
        n_ids += 1
        chunk_len += len(line) + 1
        if n_ids >= batch_size:
            jobs_for_file.append(f'<job file="{idfile}" offs="{offs}" len="{chunk_len}" subject="{subject}"></job>')
            offs += chunk_len
            chunk_len = 0
            n_ids = 0
            if len(jobs_for_file) >= batches_per_jobfile:
                with open(f"jobs/job.{job_num}", "wt") as jobfile:
                    print(os.linesep.join(jobs_for_file), file=jobfile)
                    job_num += 1
                    jobs_for_file = []
    if n_ids > 0:
        jobs_for_file.append(f'<job file="{idfile}" offs="{offs}" len="{chunk_len}" subject="{subject}"></job>')
    if jobs_for_file:
        with open(f"jobs/job.{job_num}", "wt") as jobfile:
            print(os.linesep.join(jobs_for_file), file=jobfile)
    print(batches_per_jobfile)
    """
    stub:
    """
    mkdir -p jobs
    touch jobs/1
    touch jobs/2
    touch jobs/3
    touch jobs/4
    touch ${sampleID}.ids
    touch ${sampleID}.LDS
    echo 1
    """
}


process minimap2_wnode {
    label 'long_job'
    input:
        path genome_fasta, stageAs: "genome/*"
        path genome_index // minimap2 index, addressed indirectly via job file
        path gencoll
        tuple val(sampleID), path(job), path(reads_fasta, stageAs: "reads/*"), path(ids), path(reads_LDS), val(start_job_id)
        val max_intron
        val parameters
    output:
        path "alignments/*", emit: "alignments"
    script:
        job_num = job.toString().tokenize('.').last().toInteger()
        String minimap2_params = merge_params("", parameters, "minimap2-params")
        String minimap2_wnode_params =  merge_params("-max-intron ${max_intron}", parameters, "minimap2_wnode") + ' -minimap2-params "' + minimap2_params + '"'
    """
    mkdir -p interim
    mkdir -p tmp
    mkdir -p asncache
    prime_cache -cache asncache/ -ifmt fasta -i genome/* -split-sequences
    minimap2_wnode -minimap2-executable `which minimap2` -filter-executable `which exon_selector` -start-job-id $start_job_id -input-jobs $job -nogenbank -lds2 $reads_LDS -asn-cache asncache -gc $gencoll -work-area tmp -O interim $minimap2_wnode_params
    mkdir -p alignments
    for f in interim/*; do
        if [ -f \$f ]; then
            mv \$f alignments/${sampleID}_${job_num}_\$(basename \$f)
        fi
    done
    """
    stub:
        job_num = job.toString().tokenize('.').last().toInteger()
        String minimap2_params = merge_params("", parameters, "minimap2-params")
        String minimap2_wnode_params =  merge_params("-max-intron ${max_intron}", parameters, "minimap2_wnode") + ' -minimap2-params "' + minimap2_params + '"'
        println("Effective minimap2_wnode parameters: $minimap2_wnode_params")
    """
    mkdir -p alignments
    touch alignments/${sampleID}_${job_num}.sam
    """
}


process gpx_qdump {
    input:
        path files, stageAs: "alignments/*"
    output:
        path "align.asnb", emit: "alignments"
    script:
    """
    # NB: only combination of -default-output-name and -output @.suffix
    # dumps all of jobs from the job result files
    gpx_qdump -default-output-name align -unzip '*' -input-path alignments -output @.asnb
    """
    stub:
    """
    touch align.asnb
    """
}
