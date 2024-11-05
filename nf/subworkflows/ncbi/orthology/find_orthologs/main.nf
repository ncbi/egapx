#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'

workflow find_orthologs {
    take:
        input_gencoll_asn     //path: gencoll
        ref_gencoll_asn     //path: gencoll
        input_annotations       //path: annotations
        ref_annotations        //path: annotations
        prot_hits
        blastdb
        input_proteins_asn
        ref_proteins_asn
        input_genome_asn
        ref_genome_asn
        parameters  // Map : extra parameter and parameter update
    main:
        default_params = "-check_exome -annots1_serial_type Seq-annot -annots2_serial_type Seq-annot "
        effective_params = merge_params(default_params, parameters, 'find_orthologs')
        run_find_orthologs(input_gencoll_asn, ref_gencoll_asn, input_annotations, ref_annotations, prot_hits, blastdb, 
                input_proteins_asn, ref_proteins_asn, input_genome_asn, ref_genome_asn, effective_params)

    emit:
        orthologs = run_find_orthologs.out.orthologs
        stats = run_find_orthologs.out.stats
        all = run_find_orthologs.out.all
}


process run_find_orthologs {
    input:
        path input_gencoll_asn     
        path ref_gencoll_asn, stageAs: 'input/ref_gencoll.asn'
        path input_annotations   
        path ref_annotations
        path prot_hits
        path blastdb
        path input_proteins_asn
        path ref_proteins_asn
        path input_genome_asn
        path ref_genome_asn, stageAs: 'input/ref_genome.asn'
        val parameters
    output:
        path ('output/orthologs.rpt'), emit: 'orthologs'
        path ('output/stats.xml'), emit: 'stats'
        path ('output/*'), emit: 'all'
    script:
    """
    mkdir -p output
    echo "${input_annotations.join('\n')}"  > annotations1.mft
    echo "${ref_annotations.join('\n')}"  > annotations2.mft
    str=""
    if [ -z "$blastdb" ]
    then
        mkdir -p ./asncache/
        prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i $ref_proteins_asn  -oseq-ids /dev/null -split-sequences
        prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i $input_genome_asn  -oseq-ids /dev/null -split-sequences
        prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i $input_proteins_asn  -oseq-ids /dev/null -split-sequences
        prime_cache -cache ./asncache/ -ifmt asnb-seq-entry  -i input/ref_genome.asn  -oseq-ids /dev/null -split-sequences
        str="-asn-cache ./asncache/  -prot_hits_serial_type Seq-align-set"
    else
        str="-blastdb $blastdb/prot.blastdb,$blastdb/nucl.blastdb -prot_hits_serial_type Seq-align"
    fi
    find_orthologs $parameters -gc1 $input_gencoll_asn -gc2 input/ref_gencoll.asn -annots1 annotations1.mft  -annots2 annotations2.mft  \
                 \$str  -o_orthologs output/orthologs.rpt   -prot_hits $prot_hits \
                 -o_stats output/stats.xml -nogenbank

    """
    stub:
    """
    mkdir -p output
    touch output/stats.xml
    touch output/orthologs.rpt
    """
}




//ref_prot_url='https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/ortholog_references/9606/current/GCF_000001405.40_GRCh38.p14_protein.faa.gz'
//ref_genf_url='https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/ortholog_references/9606/current/GCF_000001405.40_GRCh38.p14_genomic.fna.gz'
//ref_geng_url='https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/ortholog_references/9606/current/GCF_000001405.40_GRCh38.p14_genomic.gff.gz' 
//ref_name_url='https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/ortholog_references/9606/name_from_ortholog.rpt' 
process fetch_ortholog_references  {
    input:
        //path ortho_files      // map with file sources
        path genomic_fna
        path genomic_gff
        path protein_faa
        path name_from_rpt
    output:
        path "output/ref_protein.faa", emit: "ref_protein_faa"
        path "output/ref_genomic.fna", emit: "ref_genomic_fna"
        path "output/ref_genomic.gff.asnt.gz", emit: "annot_file"
        path "output/name_from_ortholog.rpt", emit: "name_from_ortholog"
    script:
    """
        mkdir -p output
        
        cp ${genomic_fna} output/ref_genomic.fna.gz
        cp ${protein_faa} output/ref_protein.faa.gz
        cp ${name_from_rpt} output/name_from_ortholog.rpt

        gunzip output/ref_genomic.fna.gz
        gunzip output/ref_protein.faa.gz
        zcat ${genomic_gff} | multireader -format gff3 | gzip -c > output/ref_genomic.gff.asnt.gz
    """   
    stub:
    """
        mkdir -p output
        touch output/ref_protein.faa
        touch output/ref_genomic.fna
        touch output/ref_genomic.gff.asnt.gz
        touch output/name_from_ortholog.rpt
        
        ##cp ${genomic_fna} output/ref_genomic.fna
        ##cp ${genomic_gff} output/ref_genomic.gff
        ##cp ${protein_faa} output/ref_protein.faa
        ##cp ${name_from_rpt} output/name_from_ortholog.rpt
    """
}


/*
#!/usr/bin/bash
set -exuo pipefail

gc1=GCF_000364345.1   # macaca       ([Q]uery)
gc2=GCF_000001405.40  # homo sapiens ([S]ubject)

inp_dir1=/am/ftp-genomes/all/GCF/000/364/345/GCF_000364345.1_Macaca_fascicularis_5.0/
inp_dir2=/am/ftp-genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/


cleanup() {
    kill -- -$$
}
trap cleanup SIGTERM


mkdir -p data

# Required by find_orthologs
gc_get_assembly -acc $gc1 -level 2 -o data/$gc1.gc.asn
gc_get_assembly -acc $gc2 -level 2 -o data/$gc2.gc.asn

# Convert GFF3 to ASN.1
zcat $inp_dir1/*_genomic.gff.gz | multireader -format gff3 | gzip -c > data/$gc1.annot.asnt.gz
zcat $inp_dir2/*_genomic.gff.gz | multireader -format gff3 | gzip -c > data/$gc2.annot.asnt.gz

# Genomic sequence required for aligning nucleotide neighborhood beween putative pairs of orthologs (-check_exome)                                                                                      
# Protein sequence required for computing ortholog-specific scores.
# (takes about 5 minutes)
zcat $inp_dir1/*_genomic.fna.gz $inp_dir2/*_genomic.fna.gz | makeblastdb -dbtype nucl -input_type fasta -parse_seqids -title nucs  -out data/nucl.blastdb                                               
zcat $inp_dir1/*_protein.faa.gz $inp_dir2/*_protein.faa.gz | makeblastdb -dbtype prot -input_type fasta -parse_seqids -title prots -out data/prot.blastdb                                               

# Make diamond-db for subject sequences
zcat $inp_dir2/*_protein.faa.gz | diamond makedb --db $gc2.prots

# Compute protein hits (takes about 3 minutes with 96 CPUs)
time zcat $inp_dir1/*_protein.faa.gz |
    nice -n19 diamond blastp --db ./$gc2.prots.dmnd --very-sensitive --sam-query-len --outfmt sam |
    sam2asn -diamond -align-type prot-to-prot -nogenbank |
    gzip -c > data/$gc1-$gc2.prot-hits.seq-align.asnb.gz

# Compute orthologs (takes about 10 minutes)
zcat data/$gc1-$gc2.prot-hits.seq-align.asnb.gz |
    time find_orthologs -prot_hits - -prot_hits_serial_type Seq-align \
        -gc1 data/$gc1.gc.asn \
        -gc2 data/$gc2.gc.asn \
        -annots1_serial_type Seq-annot -annots1 <(ls data/$gc1.annot.asnt.gz) \
        -annots2_serial_type Seq-annot -annots2 <(ls data/$gc2.annot.asnt.gz) \
        -check_exome \
        -blastdb data/prot.blastdb,data/nucl.blastdb \
        -o_orthologs orthologs.rpt \
        -nogenbank

*/
