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


workflow minimap2_fasta {
    take:
        genome_fasta
        genome_index // minimap2 index
        gencoll
        reads       // channel: FASTA file pairs (may or may not be pre-split)
        max_intron  // max intron length
        parameters  // Map : extra parameter and parameter update
    main:
        // Convert single-value queue channels to value channels so they can be reused
        def genome_fasta_val = genome_fasta.first()
        def genome_index_val = genome_index.first()
        def gencoll_val = gencoll.first()
        
        def split_reads = parameters.get('split', 0) as int
        
        // Optionally split files if split parameter is provided and > 1
        if (split_reads > 1) {
            split_result = split_fasta_files(reads, split_reads)
            // Emit individual files with computed job_id: part_number + (end - 1) * num_parts
            // Capture split_reads in a final variable for closure
            final int num_parts = split_reads
            split_individual = split_result.flatMap { sra, files_1, files_2 ->
                def list_1 = files_1 instanceof List ? files_1.sort() : (files_1 ? [files_1] : [])
                def list_2 = files_2 instanceof List ? files_2.sort() : (files_2 ? [files_2] : [])
                
                def result = []
                // .1. files: job_id = part_number (end=1)
                for (int idx = 0; idx < list_1.size(); idx++) {
                    result.add([sra, list_1[idx], idx + 1])
                }
                // .2. files: job_id = part_number + num_parts (end=2)
                for (int idx = 0; idx < list_2.size(); idx++) {
                    result.add([sra, list_2[idx], idx + 1 + num_parts])
                }
                return result
            }
            individual_files = split_individual
        } else {
            // Process file pairs directly - handle both (sra, files) and (sra, files_1, files_2) tuples
            // For pre-split files (from fetch_sra_fasta), extract part number from filename
            individual_files = reads.flatMap { tuple_data ->
                def sra = tuple_data[0]
                def files_1, files_2
                
                // Handle both 2-element and 3-element tuple structures
                if (tuple_data.size() >= 3) {
                    // (sra, files_1, files_2) from fetch_sra_fasta
                    files_1 = tuple_data[1] instanceof List ? tuple_data[1] : (tuple_data[1] ? [tuple_data[1]] : [])
                    files_2 = tuple_data[2] instanceof List ? tuple_data[2] : (tuple_data[2] ? [tuple_data[2]] : [])
                } else {
                    // (sra, files) - single file or list of files
                    def all_files = tuple_data[1] instanceof List ? tuple_data[1] : [tuple_data[1]]
                    files_1 = all_files.findAll { f -> f.getName() =~ /[._]1[._]/ || !(f.getName() =~ /[._]2[._]/) }
                    files_2 = all_files.findAll { f -> f.getName() =~ /[._]2[._]/ }
                }
                
                // Helper to extract part number from filename like "SRR_1.part_001.fasta.gz"
                def extractPartNumber = { fname ->
                    def matcher = fname =~ /\.part_(\d+)\./
                    if (matcher) {
                        return matcher[0][1] as int
                    }
                    return 1  // default for non-split files
                }
                
                def list_1 = files_1.sort()
                def list_2 = files_2.sort()
                def num_parts = Math.max(list_1.size(), list_2.size())
                // If single file per tuple (from fetch_sra_fasta pre-split), get num_parts from filename
                if (num_parts == 1 && list_1.size() == 1) {
                    def fname = list_1[0].getName()
                    if (fname =~ /\.part_(\d+)\./) {
                        // This is a pre-split file, use part number directly
                        def part_num = extractPartNumber(fname)
                        def job_id = fname =~ /[._]2[._]/ ? part_num + 16 : part_num  // Assume 16 parts for _2 offset
                        return [[sra, list_1[0], job_id]]
                    }
                }
                if (num_parts == 1 && list_2.size() == 1) {
                    def fname = list_2[0].getName()
                    if (fname =~ /\.part_(\d+)\./) {
                        def part_num = extractPartNumber(fname)
                        return [[sra, list_2[0], part_num + 16]]  // _2 files get offset
                    }
                }
                if (num_parts == 0) num_parts = 1
                
                def result = []
                // _1 files: job_id = part_number (end=1)
                for (int idx = 0; idx < list_1.size(); idx++) {
                    def part_num = extractPartNumber(list_1[idx].getName())
                    result.add([sra, list_1[idx], part_num])
                }
                // _2 files: job_id = part_number + num_parts (end=2)
                for (int idx = 0; idx < list_2.size(); idx++) {
                    def part_num = extractPartNumber(list_2[idx].getName())
                    result.add([sra, list_2[idx], part_num + num_parts])
                }
                return result
            }
        }
        // individual_files.view { "individual_files: $it" }
        alignments = minimap2_wnode_fasta(genome_fasta_val, genome_index_val, gencoll_val, individual_files, max_intron, parameters)
        gpx_qdump(alignments.collect())
    emit:
        alignments = gpx_qdump.out.alignments
}


process split_fasta_files {
    label 'single_cpu'
    label 'small_mem'
    input:
        tuple val(sampleID), path(fasta_rna_files, stageAs: "reads/*")
        val split_parts
    output:
        // Output patterns: .1. files (required) and .2. files (optional for paired-end)
        tuple val(sampleID), path('output/*.1.*', arity: '1..*'), path('output/*.2.*', arity: '0..*')
    script:
    """
    mkdir -p output
    
    # Rename files to use .1. and .2. convention if they use _1 and _2
    for f in reads/*; do
        [ -e "\$f" ] || continue
        case "\$f" in
            *_1.*) newname=\$(echo "\$f" | sed 's/_1\\./.1./'); mv "\$f" "\$newname" ;;
            *_2.*) newname=\$(echo "\$f" | sed 's/_2\\./.2./'); mv "\$f" "\$newname" ;;
        esac
    done
    
    # For files without .1. or .2. (single-end), rename to .1. convention
    for f in reads/*; do
        [ -e "\$f" ] || continue
        base=\$(basename "\$f")
        # Skip if already has .1. or .2. in name
        if [[ "\$base" == *.1.* ]] || [[ "\$base" == *.2.* ]]; then
            continue
        fi
        # Handle compression extensions
        comp_ext=""
        name="\$base"
        if [[ "\$name" == *.gz ]]; then
            comp_ext=".gz"
            name="\${name%.gz}"
        elif [[ "\$name" == *.zst ]]; then
            comp_ext=".zst"
            name="\${name%.zst}"
        fi
        ext="\${name##*.}"
        stem="\${name%.*}"
        newname="\${stem}.1.\${ext}\${comp_ext}"
        mv "reads/\$base" "reads/\$newname"
    done
    
    # Split each file - seqkit handles compressed input automatically
    # Use .gz extension to force gzip-compressed output (seqkit preserves format extension)
    for f in reads/*; do
        [ -e "\$f" ] || continue
        seqkit split2 -p ${split_parts} -O output --extension .gz "\$f"
    done
    """
    stub:
    """
    mkdir -p output
    touch output/${sampleID}.1.part_001.fasta.gz
    touch output/${sampleID}.2.part_001.fasta.gz
    """
}


process minimap2_wnode_fasta {
    maxForks 16
    label 'long_job'
    label 'multi_node'
    label 'large_mem'
    input:
        path genome_fasta, stageAs: "genome/*"
        path genome_index // minimap2 index
        path gencoll
        tuple val(sampleID), path(reads_fasta, stageAs: "reads/*"), val(start_job_id)
        val max_intron
        val parameters
    output:
        path "alignments/*", emit: "alignments"
    script:
        def nthreads=task.ext.threads
        String minimap2_params = merge_params("-t ${nthreads}", parameters, "minimap2-params")
        String minimap2_wnode_params =  merge_params("-max-intron ${max_intron}", parameters, "minimap2_wnode") +
                                        ' ' + merge_params("", parameters, "minimap2_wnode_fasta") +
                                        ' -minimap2-params "' + minimap2_params + '"'
    """
    # Create job.xml inline - one job entry per read file
    for f in reads/*; do
        [ -e "\$f" ] || continue
        printf '<job query="lcl|%s" subject="%s"></job>\\n' "\$f" "${genome_index}"
    done > job.xml
    
    mkdir -p tmp/asncache
    mkdir -p tmp/interim
    prime_cache -cache tmp/asncache/ -ifmt fasta -i genome/* -split-sequences
    minimap2_wnode -separate-output-by-run no -minimap2-executable `which minimap2` -filter-executable `which exon_selector` -start-job-id $start_job_id -input-jobs job.xml -nogenbank -asn-cache tmp/asncache -gc $gencoll -work-area tmp -O tmp/interim $minimap2_wnode_params
    mkdir -p alignments
    cat tmp/interim/* > alignments/minimap2_wnode.${task.index}.gpx-job.asnb
    rm -rf tmp
    """
    stub:
        String minimap2_params = merge_params("", parameters, "minimap2-params")
        String minimap2_wnode_params =  merge_params("-max-intron $max_intron", parameters, "minimap2_wnode") +
                                        ' ' + merge_params("", parameters, "minimap2_wnode_fasta") +
                                        ' -minimap2-params "' + minimap2_params + '"'
        println("Effective minimap2_wnode parameters: $minimap2_wnode_params")
    """
    mkdir -p alignments
    touch alignments/minimap2_wnode.${task.index}.gpx-job.asnb
    """
}


process gpx_qsubmit {
    label 'gpx_submitter'
    label 'small_mem'
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

    # Check for .zst files and decompress them
    reads_dir = "reads"
    for filename in os.listdir(reads_dir):
        if filename.endswith('.zst'):
            zst_path = os.path.join(reads_dir, filename)
            output_path = zst_path[:-4]  # Remove .zst extension
            
            # Skip if already decompressed
            if not os.path.exists(output_path):
                print(f"Decompressing {zst_path} -> {output_path}", flush=True)
                subprocess.run(['zstd', '-d', zst_path, '-o', output_path], check=True)
            else:
                print(f"Decompressed file already exists: {output_path}", flush=True)

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
    maxForks 50
    label 'long_job'
    label 'multi_node'
    label 'med_mem'
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
        String minimap2_params = merge_params("", parameters, "minimap2-params")
        String minimap2_wnode_params =  merge_params("-max-intron $max_intron", parameters, "minimap2_wnode")
    """
    mkdir -p tmp/asncache
    mkdir -p tmp/interim
    prime_cache -cache tmp/asncache/ -ifmt fasta -i genome/* -split-sequences

    echo "$minimap2_params" > ./minimap2.params   # see if its empty
    if [[ -s ./minimap2.params ]]; then
        minimap2_wnode -separate-output-by-run no -minimap2-executable `which minimap2` -filter-executable `which exon_selector` -start-job-id $start_job_id -input-jobs $job -nogenbank -lds2 $reads_LDS -asn-cache tmp/asncache -gc $gencoll -work-area tmp -O tmp/interim $minimap2_wnode_params -minimap2-params \"$minimap2_params\"
    else
        minimap2_wnode -separate-output-by-run no -minimap2-executable `which minimap2` -filter-executable `which exon_selector` -start-job-id $start_job_id -input-jobs $job -nogenbank -lds2 $reads_LDS -asn-cache tmp/asncache -gc $gencoll -work-area tmp -O tmp/interim $minimap2_wnode_params
    fi
    mkdir -p alignments
    cat tmp/interim/* > alignments/minimap2_wnode.${task.index}.gpx-job.asnb
    rm -rf tmp
    """
    stub:
        String minimap2_params = merge_params("", parameters, "minimap2-params")
        String minimap2_wnode_params =  merge_params("-max-intron $max_intron", parameters, "minimap2_wnode") + ' -minimap2-params "' + minimap2_params + '"'
        println("Effective minimap2_wnode parameters: $minimap2_wnode_params")
    """
    mkdir -p alignments
    touch alignments/minimap2_wnode.${task.index}.gpx-job.asnb
    """
}


process gpx_qdump {
    label 'single_cpu'
    label 'small_mem'
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
