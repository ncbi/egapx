#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow gnomon_evidence_summary {
    take:
        genome_asn
        evidence
        evidence_slices
        gn_models
        gn_models_slices
        taxid
        parameters  // Map : extra parameter and parameter update
    main:
        String gpx_qsubmit_params =  merge_params("", parameters, 'gpx_qsubmit')
        String gnomon_report_params =  merge_params("-exclusive-threshold 800000000 -max-splice-diff 5", parameters, 'gnomon_report')
        String gnomon_summary_params =  merge_params("-max-outliers 20", parameters, 'gnomon_summary')
        String gpx_qdump_params =  merge_params("-unzip '*'", parameters, 'gnomon_summary')
        //gn_models_slices|view
        scaffolds = generate_scaffolds(evidence, evidence_slices, gn_models, gn_models_slices)
        (jobs, lines_per_file) = gpx_qsubmit(scaffolds, gn_models, gpx_qsubmit_params)
        report_parts = gnomon_report(jobs.flatten(), evidence, evidence_slices, gn_models, gn_models_slices, genome_asn, lines_per_file.first(), taxid, gnomon_report_params)
        gpx_qdump(report_parts.collect(), gpx_qdump_params)
        (report, quality_report) = create_reports(gpx_qdump.out.report)
        gnomon_summary(report, quality_report,  gn_models, genome_asn, gnomon_summary_params)
    emit:
        gnomon_evidence_summary = gnomon_summary.out.gnomon_evidence_summary
        gnomon_quality_report = quality_report
        gnomon_report = report
}


process generate_scaffolds {
    input:
        path evidence
        path evidence_slices
        path gn_models
        path gn_models_slices
    output:
        path './scaffold.list', emit: 'scaffold_list'
    script:
    """
#!/usr/bin/env python3
import os
from pathlib import Path

def main():
    files = []
    lines = []
    slice_map = {}
    files = "$evidence_slices".split()
    files.extend("$gn_models_slices".split())
    scaffold = None
    for filename in files:
        filename = Path(filename.strip()).name
        input_type = -1
        scaffold = None
        #if filename.startswith("gnomon_wnode"):
        if filename.find("gnomon_wnode.out") >= 0:
            input_type = 0
        elif filename.find("evidence") >= 0:
            input_type = 1
        else:
            continue
        with open(filename, 'r') as file:
            lines = file.readlines()
        filename = lines[0].strip()
        filename = Path(filename).name
        for line in lines[1:]:
            separator = line.find('\\t')
            pos = int(line[separator+1:])
            if scaffold is not None:  # prior row
                slice_map[scaffold][input_type][1] = pos-1
            scaffold = line[:separator]
            if scaffold not in slice_map.keys():
                slice_map[scaffold] = {0:[0,0], 1:[0,0]}
            slice_map[scaffold][input_type][0] = pos

        filesize = os.path.getsize(filename)
        slice_map[scaffold][input_type][1] = filesize-1
    pairs = []
    for scaffold, values in slice_map.items():
        def get_length(a,b):
            if a == b and a == 0:
                return 0
            else:
                return (b - a)+1
        pairs.append([scaffold, get_length(values[0][0], values[0][1]) + get_length(values[1][0], values[1][1]), values])
    pairs.sort(key = lambda x : x[1], reverse = True)
    with open('scaffold.list', 'w') as file:
        for pair in pairs:
            ##file.write(pair[0]+'\\t'+str(pair[1])+'\\t'+str(pair[2])+'\\n')
            file.write(pair[0]+'\\n')

if __name__=="__main__":
    main()
    """
    
    stub:
    """
    touch scaffold.list
    """
}

process create_reports {
    input:
        path reports
    output:
        path './gnomon_report.txt', emit: 'gnomon_report'
        path './gnomon_quality_report.txt', emit: 'gnomon_quality_report'
    script:

"""
#!/usr/bin/env python3
import os

def main():

    f_quality = open('gnomon_quality_report.txt', 'w')
    f_support = open('gnomon_report.txt', 'w')
    f_support.write('#Gnomon model	Scaffold id	Evidence id	Evidence taxid	Evidence origin	Evidence PIG id	Alignment Percent Identity	Base Coverage Percentage	CDS Base Coverage Percentage	Precise splice-site support	Approximate splice-site support	Core Support	In Minimal Full Introns Support\\n')
    f_quality.write('#Gnomon model	Scaffold id	Minimal Full Support	Minimal Same-species Full Support	Minimal Full Intron Support	Minimal Same-species Full Intron Support	Average Base Same-Species Support	Smallest Base Same-Species Support	Average Intron Same-Species Support	Smallest Intron Same-Species Support	Number Introns Same-Species Supported	Ab Initio Percentage	SRS Base Support Percentage	Full intron support SRS count	Partial intron support SRS count	Non-consensus introns	Overlapping strong-signal uORFs	Non-overlapping strong-signal uORFs	SRS Base Support Percentage Unambiguous\\n')
    file = None
    with open('$reports', 'r') as file:
        for line in file:
            line = line.strip()
            if line == 'QUALITY':
                file = f_quality
                continue
            if line == 'SUPPORT':
                file = f_support
                continue
            if file is None:
                print("file not well formatted")
                exit()
            file.write(line + '\\n')

if __name__=="__main__":
    main()
"""
    stub:
    """
    touch gnomon_quality_report.txt
    touch gnomon_report.txt
    """
}



process gpx_qsubmit {
    input:
        path scaffolds
        path gn_models
        val params
    output:
        path "job.*"
        env lines_per_file
    script:
        njobs=16
    """
    echo $scaffolds > scaffolds.mft
    mkdir -p ./asncache/
    prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ${gn_models} -oseq-ids spids -split-sequences

    gpx_qsubmit $params -ids-manifest scaffolds.mft -o jobs -nogenbank -asn-cache ./asncache/  -keep-input-order
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
    stub:
        njobs=16
    """
    for i in {1..$njobs}; do
        echo j.\${i} >> jobs
    done
    split -nr/$njobs jobs job. -da 3
    lines_per_file=10
    """
}

process gnomon_report {
    input:
        path jobs
        path evidence
        path evidence_slices
        path gn_models
        path gn_models_slices
        path genome_asn
        val lines_per_file
        val taxid
        val params
    output:
        path "output/*", emit: "output"
    script:
        job_num = jobs.toString().tokenize('.').last().toInteger()
    """
    njobs=`wc -l <$jobs`
    if [ \$njobs -lt 16 ]; then
        threads=\$njobs
    else
        threads=16
    fi  
    mkdir -p ./asncache/
    prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ${gn_models} -oseq-ids spids2 -split-sequences
    prime_cache -cache ./asncache/ -ifmt asn-seq-entry  -i ${genome_asn} -oseq-ids spids -split-sequences

    filename=\$(basename -- "$jobs")
    extension="\${filename##*.}"
    (( start_job_id = ((10#\$extension) * $lines_per_file) + 1 ))

    # When running multiple jobs on the cluster there is a chance that
    # several jobs will run on the same node and thus generate files
    # with the same filename. We need to avoid that to be able to stage
    # the output files for gnomon_report. We add the job file numeric
    # extension as a prefix to the filename.
    echo "${evidence_slices.join('\n')}" > input_slices.mft
    echo "${gn_models_slices.join('\n')}" >> input_slices.mft
    mkdir -p interim
    #gnomon_report $params -nogenbank -input-slices input_slices.mft -same-species $taxid -asn-cache ./asncache/  -start-job-id \$start_job_id -workers \$threads -input-jobs $jobs -O interim
    gnomon_report $params -nogenbank -input-slices input_slices.mft -asn-cache ./asncache/  -start-job-id \$start_job_id -workers \$threads -input-jobs $jobs -O interim || true
    mkdir -p output
    for f in interim/*; do
        if [ -f \$f ]; then
            mv \$f output/\${extension}_\$(basename \$f)
        fi
    done
    """
    stub:
        job_num = jobs.toString().tokenize('.').last().toInteger()
    """
    mkdir -p output
    touch output/gnomon_report.${job_num}.txt
    """
}


process gnomon_summary {
    input:
        path gnomon_report
        path gnomon_quality_report
        path gn_models
        path genome_asn
        val params
    output:
        path "All.gnomon_evidence_*", emit: "gnomon_evidence_summary"
    script:
    """
    
    # mkdir -p ./asncache/
    # prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i ${gn_models} -oseq-ids spids2 -split-sequences
    # prime_cache -cache ./asncache/ -ifmt asn-seq-entry  -i ${genome_asn} -oseq-ids spids -split-sequences
    # echo "${gn_models.join('\n')}" > models.mft
    # gnomon_summary -models models.mft -nogenbank $params  -input $gnomon_report -quality $gnomon_quality_report -asn-cache ./asncache/  -output All.gnomon_evidence_@.txt

    touch All.gnomon_evidence_bag.txt
    touch All.gnomon_evidence_pair.txt
    touch All.gnomon_evidence_stats.txt
    touch All.gnomon_evidence_taxid.txt
    """
    stub:
    """
    touch All.gnomon_evidence_bag.txt
    touch All.gnomon_evidence_pair.txt
    touch All.gnomon_evidence_stats.txt
    touch All.gnomon_evidence_taxid.txt
    """
}


process gpx_qdump {
    input:
        path files, stageAs: "inputs/*"
        val params
    output:
        path "reports.txt", emit: "report"
    script:
    """
    gpx_qdump $params -input-path inputs -o reports.txt 
    touch reports.txt 
    """
    stub:
    """
    touch reports.txt 
    """
}
