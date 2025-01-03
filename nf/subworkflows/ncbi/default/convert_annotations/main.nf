#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow convert_annotations {
    take:
        asn_file    // Channel: ASN file with annotations
        parameters  // Map : extra parameter and parameter update
                    // Use "annotwriter_gff" for GFF3 and "annotwriter_gtf" for GTF in section 'convert_annotations'
    main:
        run_converter(asn_file, merge_params("", parameters, 'annotwriter_gff'), merge_params("", parameters, 'annotwriter_gtf'))
    emit:
        genomic_gff = run_converter.out.genomic_gff
        genomic_gtf = run_converter.out.genomic_gtf
        genomic_fasta = run_converter.out.genomic_fasta
        transcripts_fasta = run_converter.out.transcripts_fasta
        cds_fasta = run_converter.out.cds_fasta
        proteins_fasta = run_converter.out.proteins_fasta
}


process run_converter {
    input:
        path asn_files, stageAs: 'asn_inputs/*'
        val  gff_params
        val  gtf_params
    output:
        path 'output/*.genomic.gff', emit: 'genomic_gff'
        path 'output/*.genomic.gtf', emit: 'genomic_gtf'
        path 'output/*.genomic.fna', emit: 'genomic_fasta'
        path 'output/*.transcripts.fna', emit: 'transcripts_fasta'
        path 'output/*.cds.fna', emit: 'cds_fasta'
        path 'output/*.proteins.faa', emit: 'proteins_fasta'
    script:
    """
    mkdir -p output
    mkdir -p tmpout
    found_afbs=(0)
    for af in asn_inputs/*
    do
        afb=\$(basename \$af)
        found_afbs+=(\${afb})
        annotwriter ${gff_params} -nogenbank -i \${af} -format gff3 -o tmpout/\${afb}.genomic.gff
        annotwriter ${gtf_params} -nogenbank -i \${af} -format gtf -o tmpout/\${afb}.genomic.gtf
        asn2fasta -nogenbank -i \${af} -nucs-only |sed -e 's/^>lcl|\\(.*\\)/>\\1/' > tmpout/\${afb}.genomic.fna
        asn2fasta -nogenbank -i \${af} -feats rna_fasta -o tmpout/\${afb}.transcripts.fna
        asn2fasta -nogenbank -i \${af} -feats fasta_cds_na -o tmpout/\${afb}.cds.fna
        asn2fasta -nogenbank -i \${af} -prots-only -o tmpout/\${afb}.proteins.faa
    done
    ##echo 'D: ' \${found_afbs[@]}
    cat `find tmpout -name g*.gff -o -name all_unannot*.genomic.gff` > output/complete.genomic.gff
    cat `find tmpout -name g*.gtf -o -name all_unannot*.genomic.gtf` > output/complete.genomic.gtf
    cat `find tmpout -name g*.genomic.fna -o -name all_unannot*.genomic.fna` > output/complete.genomic.fna
    cat `find tmpout -name g*.transcripts.fna -o -name all_unannot*.transcripts.fna` > output/complete.transcripts.fna
    cat `find tmpout -name g*.cds.fna -o -name all_unannot*.cds.fna` > output/complete.cds.fna
    cat `find tmpout -name g*.proteins.faa -o -name all_unannot*.proteins.faa` > output/complete.proteins.faa
    rm tmpout/*
    touch output/complete.genomic.gff
    touch output/complete.genomic.gtf
    touch output/complete.genomic.fna
    touch output/complete.transcripts.fna
    touch output/complete.cds.fna
    touch output/complete.proteins.faa
    """

    stub:
    """
    mkdir -p output
    echo "Genomic GFF"      > output/complete.genomic.gff
    echo "Genomic GTF"      > output/complete.genomic.gtf
    echo "Genomic FASTA"    > output/complete.genomic.fna
    echo "Transcript FASTA" > output/complete.transcripts.fna
    echo "CDS FASTA"        > output/complete.cds.fna
    echo "Protein FASTA"    > output/complete.proteins.faa
    """
}
