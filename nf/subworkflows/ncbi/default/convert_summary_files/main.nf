#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow convert_summary_files {
    take:
        report
        quality_report
        locustypes
        locus_tag_prefix
    main:
        (report, quality_report) = run_convert_summary_files(report, quality_report,  locustypes, locus_tag_prefix)
    emit:
        gnomon_quality_report = quality_report
        gnomon_report = report
}


process run_convert_summary_files {
    input:
        path report
        path quality_report
        path locustypes
        val locus_tag_prefix 
    output:
        path "./new.$report", emit: 'report'
        path "./new.$quality_report", emit: 'quality_report'       
    script:
    """
#!/usr/bin/env python3
import os
fd = open('$locustypes')
lines = fd.read().splitlines()
seq_id_to_transcriptID = {}

locus_tag_prefix = '$locus_tag_prefix' or 'egapxtmp'
    
for line in lines[1:]: #skip the first line
    line = line.split('\\t')

    if len(line) < 15:
        continue

    feat_type = line[3]
    pre_product = line[5]

    if feat_type[-3:] != 'RNA' or not pre_product:  # e.g. tRNAs
        continue

    sep_index = pre_product[pre_product.index('|') + 1:].index('|') + pre_product.index('|') + 2
    seq_id = pre_product[sep_index:-2]

    assigned_product = line[7]

    if assigned_product[0:3] not in ("XM_", "XR_"):
        continue

    geneID = line[4]

    if geneID == '-1':
        continue

    if int(geneID) > 990000000:
        geneID = geneID[3:]

    locus_tag = locus_tag_prefix + '_' + geneID
    assigned_name = line[14]
    trans = 'transcript variant X'
    pos = assigned_name.find(trans)

    #if (pos := assigned_name.find(trans)) == -1:

    if pos == -1:
        transcript_id = locus_tag + "-R1"
    else:
        transcript_id = locus_tag + "-R" + assigned_name[pos + len(trans) :]

    seq_id_to_transcriptID[seq_id] = transcript_id


def replace_identifier(seq_id_to_transcriptID, file):
    newfile = "new." + file
    is_header = True

    with open(file) as infile, open(newfile, 'w') as outfile:
        for line in infile:
            line = line.split('\\t')

            if is_header:
                line[0] = line[0][1:]
                line.insert(0,'#transcript_id ')
                line = "\\t".join(line)
                outfile.write(line)
                is_header = False
                continue

            if len(line) < 5:
                continue

            seq_id = line[0]
            transcript_id = seq_id_to_transcriptID.get(seq_id, '')
            line.insert(0,transcript_id)
            line = "\\t".join(line)
            outfile.write(line)

replace_identifier(seq_id_to_transcriptID, "$report")
replace_identifier(seq_id_to_transcriptID, "$quality_report")
    """
    
    stub:
    """
    touch new.$report
    touch new.$quality_report
    """
}
