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
        //def basename = asn_file.baseName.toString()
        def basename = asn_files.first().baseName.toString()
    """
    echo "${asn_files.join('\n')}" > ${basename}.mft
    mkdir -p output
    ##if [ -s ${asn_files} ]; then
        mkdir -p tmpout
        for af in ${asn_files}
        do
          afb=\$(basename \$af)
          annotwriter ${gff_params} -nogenbank -i \${af} -format gff3 -o tmpout/\${afb}.genomic.gff
          annotwriter ${gtf_params} -nogenbank -i \${af} -format gtf -o tmpout/\${afb}.genomic.gtf
          asn2fasta -nogenbank -i \${af} -nucs-only |sed -e 's/^>lcl|\\(.*\\)/>\\1/' > tmpout/\${afb}.genomic.fna
          asn2fasta -nogenbank -i \${af} -feats rna_fasta -o tmpout/\${afb}.transcripts.fna
          asn2fasta -nogenbank -i \${af} -feats fasta_cds_na -o tmpout/\${afb}.cds.fna
          asn2fasta -nogenbank -i \${af} -prots-only -o tmpout/\${afb}.proteins.faa
        done
        cat tmpout/*.gff > output/complete.genomic.gff
        cat tmpout/*.gtf > output/complete.genomic.gtf
        cat tmpout/*.genomic.fna > output/complete.genomic.fna
        cat tmpout/*.transcripts.fna > output/complete.transcripts.fna
        cat tmpout/*.cds.fna > output/complete.cds.fna
        cat tmpout/*.proteins.faa > output/complete.proteins.faa
        rm tmpout/*
      
        ##annotwriter ${gff_params} -nogenbank -i ${asn_files} -format gff3 -o output/${basename}.genomic.gff
        ##annotwriter ${gtf_params} -nogenbank -i ${asn_files} -format gtf -o output/${basename}.genomic.gtf
        ##asn2fasta -nogenbank -nucs-only -indir asn_inputs  -o - |sed -e 's/^>lcl|\\(.*\\)/>\\1/' >output/${basename}.genomic.fna
        ##asn2fasta -nogenbank -feats rna_fasta        -indir asn_inputs -o output/${basename}.transcripts.fna
        ##asn2fasta -nogenbank -feats fasta_cds_na -i  -indir asn_inputs -o output/${basename}.cds.fna
        ##asn2fasta -nogenbank -prots-only -i -indir asn_inputs -o output/${basename}.proteins.faa
    ##else
    ##    touch output/${basename}.genomic.gff
    ##    touch output/${basename}.genomic.gtf
    ##    touch output/${basename}.genomic.fna
    ##    touch output/${basename}.transcripts.fna
    ##    touch output/${basename}.cds.fna
    ##    touch output/${basename}.proteins.faa
    ##fi
    """

    stub:
        def basename = asn_files.first().baseName.toString()
        print(asn_files)
        print(basename)
    """
    mkdir -p output
    echo "Genomic GFF"    > output/${basename}.genomic.gff
    echo "Genomic GTF"    > output/${basename}.genomic.gtf
    echo "Genomic FASTA"  > output/${basename}.genomic.fna
    echo "Transcript FASTA" > output/${basename}.transcripts.fna
    echo "CDS FASTA" > output/${basename}.cds.fna
    echo "Protein FASTA"   > output/${basename}.proteins.faa
    """
}
