#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'

///final_asn 
//  -scaffolds final_asn_markup.8279392/inp/scaffold.mft 
//  -chromosomes final_asn_markup.8279392/inp/chromosome.mft 
//  -annots final_asn_markup.8279392/inp/annotation.mft 
//  -gene_weights final_asn_markup.8279392/inp/gene_weight.mft 
//  -locus_lnk final_asn_markup.8279392/inp/gene_assignment.mft 
//  -out_dir final_asn_markup.8279392/out 
//  -lsq_dir final_asn_markup.8279392/var/idx 
//  -gencoll-asn passthrough_gcaccess.8279042/out/gencoll.asn 
//  -asn-cache sequence_cache 
//  -nogenbank

// the scaf and chrom files are all in bare_scaffold_asn format, individual chroms, scafs bundled in per Asm and AsmUnit files
// final_asn task also runs asn_cleanup and asnvalidate on everything, examples from runlog below
// and it calls  asn_stats, asnval2gbitem (?) , gbproject, and load_finl_asn_tracking_data after it all, some of which will probably be discarded by egapx

///prime_cache -cache final_asn_markup.8279392/var/nucprots.cache -input-manifest final_asn_markup.8279392/out/nucprot.mft -ifmt asn-seq-entry
///asn_cleanup -basic -i final_asn_markup.8279392/out/raw/scaf/GCF_030936175.1/all_unannotated.asn -o final_asn_markup.8279392/out/scaf/GCF_030936175.1/all_unannotated.asn
///asn_cleanup -basic -i final_asn_markup.8279392/out/raw/scaf/GCF_030936175.1/asm_Primary_Assembly_1.cat.asn -o final_asn_markup.8279392/out/scaf/GCF_030936175.1/asm_Primary_Assembly_1.cat.asn
///asn_cleanup -basic -i final_asn_markup.8279392/out/raw/chrom/GCF_030936175.1/Chr_SUPER_1.asn -o final_asn_markup.8279392/out/chrom/GCF_030936175.1/Chr_SUPER_1.asn
///asn_cleanup -basic -i final_asn_markup.8279392/out/raw/chrom/GCF_030936175.1/Chr_SUPER_10.asn -o final_asn_markup.8279392/out/chrom/GCF_030936175.1/Chr_SUPER_10.asn

///asnvalidate -Q 0 -asn-cache final_asn_markup.8279392/var/nucprots.cache,sequence_cache -v 4 -A -X -Z -o final_asn_markup.8279392/out/val/GCF_030936175.1/asm_Primary_Assembly_1.cat.val -i final_asn_markup.8279392/out/scaf/GCF_030936175.1/asm_Primary_Assembly_1.cat.asn
///asnvalidate -Q 0 -asn-cache final_asn_markup.8279392/var/nucprots.cache,sequence_cache -v 4 -A -X -Z -o final_asn_markup.8279392/out/val/GCF_030936175.1/asm_Primary_Assembly_1001.cat.val -i final_asn_markup.8279392/out/scaf/GCF_030936175.1/asm_Primary_Assembly_1001.cat.asn
///asnvalidate -Q 0 -asn-cache final_asn_markup.8279392/var/nucprots.cache,sequence_cache -v 4 -A -X -Z -o final_asn_markup.8279392/out/val/GCF_030936175.1/Chr_SUPER_1.val -i final_asn_markup.8279392/out/chrom/GCF_030936175.1/Chr_SUPER_1.asn
///asnvalidate -Q 0 -asn-cache sequence_cache -v 4 -A -X -Z -o final_asn_markup.8279392/out/val/GCF_030936175.1/all_nucprots.scaf.val -i final_asn_markup.8279392/out/scaf/GCF_030936175.1/all_nucprots.asn
///asnvalidate -Q 0 -asn-cache sequence_cache -v 4 -A -X -Z -o final_asn_markup.8279392/out/val/GCF_030936175.1/all_nucprots.chrom.val -i final_asn_markup.8279392/out/chrom/GCF_030936175.1/all_nucprots.asn

///asn_stats -input-manifest final_asn_markup.8279392/var/joint.mft -nucprot-manifest final_asn_markup.8279392/out/nucprot.mft -o final_asn_markup.8279392/out/feature_counts.txt -counts-xml-output final_asn_markup.8279392/out/feature_counts.xml -stats-xml-output final_asn_markup.8279392/out/feature_stats.xml -t -break-by assembly-unit -asn-cache sequence_cache -gencoll-asn passthrough_gcaccess.8279042/out/gencoll.asn
///asnval2gbitem -t -asn-cache sequence_cache -asnval-path final_asn_markup.8279392/out/val/GCF_030936175.1 -scaffold-manifest final_asn_markup.8279392/out/GCF_030936175.1.scaffolds.mft -chromosome-manifest final_asn_markup.8279392/out/GCF_030936175.1.chromosomes.mft -nucprot-manifest final_asn_markup.8279392/out/GCF_030936175.1.nucprots.mft -o final_asn_markup.8279392/var/gbench/GCF_030936175.1
///gbproject -collapse final_asn_markup.8279392/var/gbench/GCF_030936175.1 -o final_asn_markup.8279392/out/GCF_030936175.1.gbp
///load_final_asn_tracking_data -feature-counts-xml final_asn_markup.8279392/out/feature_counts.xml -length-stats-xml final_asn_markup.8279392/out/feature_stats.xml -validation-xml final_asn_markup.8279392/out/annot.val.xml


workflow final_asn_markup {
    take:
        gencoll_asn
        genome_asn
        scaffolds // asn seqentrs
        chromosomes // asn seqentrysseqids
        annots // asnt seq-annots
        locus_link // rpt from locus_link
        locustypes // tsv from locus_link
        parameters  // Map : extra parameter and parameter update
    main:
        params = merge_params("", parameters, 'final_asn')

        final_asn(gencoll_asn, genome_asn, scaffolds, chromosomes, annots, locus_link, locustypes, params)

    emit:
        outputs = final_asn.out.all
        to_convert = final_asn.out.to_convert
        validated = final_asn.out.validated
        stats = final_asn.out.stats
        annotated_genome_asn = final_asn.out.annotated_genome_asn
        annotation_data_comment = final_asn.out.annotation_data_comment
}


process final_asn {
    input:
        path gencoll_asn, stageAs: 'gencoll.asn'
        path genome_asn, stageAs: 'genome/*'
        path scaffolds, stageAs: 'scaffolds' // asn seqentry
        path chromosomes, stageAs: 'chromosomes' // asn seqentry
        path annots, stageAs: 'annots/*'  // asnt seq-annots
        path locus_link // tsv rpt
        path locustypes // tsv
        val params
    output:
        path "output/*", emit: "all"
        path "output/scaf/EGAPx_Test_Assembly/*.asn", emit: "to_convert"
        path "output/val/EGAPx_Test_Assembly/*", emit: "validated"
        path "output/stats/*", emit: "stats"
        path "output/annotated_genome.asn", emit: "annotated_genome_asn"
        path "output/annotation_data.cmt", emit: "annotation_data_comment"
    script:
    """
    mkdir -p output
    mkdir -p asncache
    mkdir -p 'EGAPx_Test_Assembly'

    prime_cache -cache ./asncache/ -ifmt asn-seq-entry  -i $genome_asn  -oseq-ids cached_ids  -split-sequences
    concat_seqentries -cache ./asncache/ -o "./EGAPx_Test_Assembly/genome.asnb.gz"
    asn_translator -gzip -i "./EGAPx_Test_Assembly/genome.asnb.gz"  -o "./EGAPx_Test_Assembly/genome.asnt" 

    echo "./EGAPx_Test_Assembly/genome.asnt" > ./scaffold.mft
    touch ./chromosome.mft
    ls -1 annots/* > ./annots.mft
    echo $locus_link > ./locus_link.mft
    echo $locustypes > ./locus_types.mft
    
    echo "" > ./gene_weights.mft

    ##lds2_indexer -source genome/ -db LDS2
    ## prime_cache
    # EXCEPTION_STACK_TRACE_LEVEL=Warning DEBUG_STACK_TRACE_LEVEL=Warning DIAG_POST_LEVEL=Trace

    final_asn $params -egapx -nogenbank  -gencoll-asn $gencoll_asn -asn-cache ./asncache/  \
        -scaffolds ./scaffold.mft  -chromosomes ./chromosome.mft  \
        -gene_weights ./gene_weights.mft  \
        -annots ./annots.mft -locus_lnk ./locus_link.mft -locus_types ./locus_types.mft \
        -S NONE -genbank-mode -out_dir  ./output/

    mkdir -p raw/scaf
    mv ./output/scaf/EGAPx_Test_Assembly/*.asn ./raw/scaf
    for f in ./raw/scaf/*.asn; do
        of=./output/scaf/EGAPx_Test_Assembly/`basename \$f`
        asn_cleanup -basic -i \$f -o \$of
        cat \$of >> output/annotated_genome.asn
    done

    # NB if (when) chromosomes is not empty the same logic should be applied to chrom directroies
    if [ -s ./output/chrom/EGAPx_Test_Assembly/*.asn ]; then
        mkdir -p raw/chrom
        mv ./output/chrom/EGAPx_Test_Assembly/*.asn ./raw/chrom
        for f in ./raw/chrom/*.asn; do
            of=./output/chrom/EGAPx_Test_Assembly/`basename \$f`
            asn_cleanup -basic -i \$f -o \$of
            cat \$of >> output/annotated_genome.asn
        done
    fi

    mkdir -p output/val/EGAPx_Test_Assembly
    for f in ./output/scaf/EGAPx_Test_Assembly/*.asn; do
        asnvalidate -Q 0 -asn-cache ./asncache/ -v 4 -A -X -Z -o ./output/val/EGAPx_Test_Assembly/`basename \$f .asn`.val -i \$f
    done

    # joint manifest is scaffolds, chromosomes, and organelles (not implemented here)
    # take it from annotated_genome.asn
    echo "./output/annotated_genome.asn" > ./joint.mft

    mkdir -p output/stats
    asn_stats -input-manifest ./joint.mft -o output/stats/feature_counts.txt -counts-xml-output output/stats/feature_counts.xml -stats-xml-output output/stats/feature_stats.xml -t -break-by assembly-unit -asn-cache ./asncache/ -gencoll-asn $gencoll_asn -genbank-mode
    """
    stub:
    """
    mkdir -p output/ACCEPT
    echo "1" > output/ACCEPT/something.asn
    mkdir -p output/scaf/EGAPx_Test_Assembly/
    echo "1" > output/scaf/EGAPx_Test_Assembly/genome.asn
    mkdir -p output/val/EGAPx_Test_Assembly/
    echo "1" > output/val/EGAPx_Test_Assembly/genome.val
    mkdir -p output/stats
    echo "1" > output/stats/feature_counts.txt
    echo "1" > output/annotated_genome.asn
    echo "1" > output/annotation_data.cmt


    echo "1" > output/final_asn.log
    """
}
