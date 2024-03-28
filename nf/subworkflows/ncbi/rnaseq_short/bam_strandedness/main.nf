#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow bam_strandedness {
    take:
        ch_bam              // channel: BAM
        ch_index            // channel: BAM index
        sra_metadata_list   // list of path: files with sra metadata
        parameters          // Map : extra parameter and parameter update
    main:
        // Not yet used but supposed to be used by rnaseq_divide_by_strandedness 
        rnaseq_divide_by_strandedness_params = merge_params("-min-aligned 1000000 -min-unambiguous 200 -min-unambiguous-pct 2 -max-unambiguous-pct 100 -percentage-threshold 98", parameters, 'rnaseq_divide_by_strandedness')
        out = exec(ch_bam, ch_index, sra_metadata_list)
        partials = out.collect()
        merge(partials, rnaseq_divide_by_strandedness_params)
    emit:
        strandedness = merge.out.strandedness
        stranded_runs = merge.out.stranded_runs
        unstranded_runs = merge.out.unstranded_runs
}


process exec {
    input:
        path sorted_bam
        path sorted_bam_index
        path sra_metadata_list
    output:
        path "output/*"
    script:
        output = "output"
    script:
        // println("Checking ${sorted_bam}")
        def (_, assembly, run_part) = (sorted_bam =~ /^([^-]+)-(.+)-[^-]+$/)[0]
    """
    #!/usr/bin/env bash

    mkdir -p $output
    echo "bam file : ${sorted_bam}, index file : ${sorted_bam_index}"
    s=\$(basename ${sorted_bam})
    # remove -Aligned.out.Sorted.bam from the string
    s=\${s%-Aligned.out.Sorted.bam}
    # Split the string by '-' to get the assembly name
    IFS=- read -a parts <<< "\$s"
    assembly=\${parts[0]}
    run=\${parts[1]}

    echo "Assembly \$assembly, run \$run"
    for mate in 1 2; do
        f=\$((mate*64))
        total_count=\$(samtools view --threads 4 -f \$f -c ${sorted_bam})
        sense_count=\$(samtools view --threads 4 -f \$f -e '[XS]=="+" && flag.reverse==0' -c ${sorted_bam})
        antisense_count=\$(samtools view --threads 4 -f \$f -e '[XS]=="+" && flag.reverse==16' -c ${sorted_bam})
        printf "\$run\\t\$mate\\t\$total_count\\t\$sense_count\\t\$antisense_count\\n" >> $output/${run_part}.\${mate}.strandedness
    done
    """
}


process merge {
    input:
        path partials
        val params
    output:
        path "run.strandedness", emit: 'strandedness'
        path "stranded.list", emit: 'stranded_runs'
        path "unstranded.list", emit: 'unstranded_runs'
    script:
    """
    #!/usr/bin/env python3
    from collections import defaultdict

    # Parse parameters:
    param_dict = {}
    k = ''
    for v in "${params}".split():
        if v[0] == "-":
            k = v[1:]
        elif k:
            param_dict[k] = int(v)
            k = ''
    percentage_threshold = param_dict["percentage-threshold"]
    min_aligned = param_dict["min-aligned"]
    min_unambiguous = param_dict["min-unambiguous"]
    min_unambiguous_pct = param_dict["min-unambiguous-pct"]
    max_unambiguous_pct = param_dict["max-unambiguous-pct"]
    

    runs_total = defaultdict(int)
    runs_sense = defaultdict(int)
    runs_antisense = defaultdict(int)
    for file in "${partials}".split():
        with open(file, "rt") as f:
            line = f.readline().strip()
            run, mate, total_count, sense_count, antisense_count= line.split('\\t')
            runs_total[(run, mate)] += int(total_count)
            runs_sense[(run, mate)] += int(sense_count)
            runs_antisense[(run, mate)] += int(antisense_count)
    strandedness_map = {}
    with open("run.strandedness", "wt") as f_strand:
        f_strand.write("#Run\\tmate\\ttotal_oount\\tsense_count\\tantisense_count\\tstrandedness\\n")
        for (run, mate) in runs_total:
            total_count = runs_total[(run, mate)]
            sense_count = runs_sense[(run, mate)]
            antisense_count = runs_antisense[(run, mate)]
            strandedness = "unstranded"

            total_unambiguous = sense_count + antisense_count
            if total_count == 0:
                strandedness="no alignments"
            elif total_count < min_aligned:
                strandedness="can't determine"
            elif total_unambiguous * 100 < total_count * min_unambiguous_pct or total_unambiguous * 100 > total_count * max_unambiguous_pct:
                strandedness="suspect"
            elif total_unambiguous < min_unambiguous:
                strandedness="can't determine"
            elif sense_count > antisense_count:
                if sense_count * 100 > total_unambiguous * percentage_threshold:
                    strandedness="sense"
                else:
                    strandedness="unstranded"
            else:
                if antisense_count * 100 > total_unambiguous * percentage_threshold:
                    strandedness="antisense"
                else:
                    strandedness="unstranded"
            if strandedness in ["sense", "antisense", "unstranded"]:
                strandedness_map[(run, mate)] = strandedness
            f_strand.write(f"{run}\\t{mate}\\t{total_count}\\t{sense_count}\\t{antisense_count}\\t{strandedness}\\n")
    stranded_set = set()
    unstranded_set = set()
    for ((run, mate), v) in strandedness_map.items():
        mate = str(3 - int(mate))
        if v == "sense" and strandedness_map[(run, mate)] == "antisense":
            stranded_set.add(run)
        elif v == "antisense" and strandedness_map[(run, mate)] == "sense":
            stranded_set.add(run)
        elif v == "unstranded" and strandedness_map[(run, mate)] == "unstranded":
            unstranded_set.add(run)

    with open("stranded.list", "wt") as f:
        for run in stranded_set:
            f.write(f"{run}\\n")
    with open("unstranded.list", "wt") as f:
        for run in unstranded_set:
            f.write(f"{run}\\n")
    """
}
