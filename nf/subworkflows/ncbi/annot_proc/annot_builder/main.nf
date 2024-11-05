#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'

// annot_builder_main collection_dir ${output}/COLLECTION accept_dir ${output}/ACCEPT conflict_dir ${output}/CONFLICT report_dir ${output}/REPORT test_dir ${output}/TEST i
//                    reftrack_attrs_manifest ${input.reftrack_attrs} loss_pct_ccds 0.0 loss_pct_refseq_ref_primary 1.0 loss_pct_refseq_alt_ref_loci 100.0 loss_pct_refseq_patches 100.0 loss_pct_refseq_other 1.0
// annot_builder_input name gnomon desc Gnomon aliases Gnomon|Chainer|PartAbInitio|FullAbInitio|Chainer_GapFilled|PartAbInitio_GapFilled 
//                     is_primary 1 input_manifest ${input.gnomon_models} model_maker gnomon2model use_secondary_support 1 keep_top_N_models 50 
//                     drop_alt_brs_overlap 1 merge_variants 1 enable_AR0050_AR0048 1 max_pct_ab_initio 50
// annot_builder -accept-output both -asn-cache ${GP_cache_dir} -conffile ${conffile} -gc-assembly-manifest ${input.gencoll_asn} -logfile ${logfile}

// not implimented, future examples
//load_annot_builder_tracking_data -stats-xml ${output}/REPORT/stats.xml -taskrun ${taskrun.id}
//annot_builder_input name bestrs desc BestRefSeq is_primary 1 input_manifest ${input.best_rs_seqalign} model_maker splign2model user_filter lxr_data.is_refseq=1
//annot_builder_input name ng desc Curated Genomic is_primary 1 input_manifest ${input.best_ng_seqalign} model_maker ng2model score_filter rank=1 exclude_subtypes CloneRef,misc_difference,STS,tRNA,variation,VariationRef exclude_types Biosrc,Pub
//annot_builder_input name imgt desc IMGT is_primary 1 input_manifest ${input.imgt} model_maker imgt2model use_secondary_support 0
//annot_builder_input name rfam desc cmsearch aliases Rfam is_primary 1 input_manifest ${input.rfam} model_maker gnomon2model
//annot_builder_input name trna desc tRNAscan-SE is_primary 1 input_manifest ${input.trna_annot} model_maker passthru
//annot_builder_input name blessed desc SelectedGeneRepresentative is_primary 0 input_manifest ${input.best_rs_seqalign} model_maker splign2model score_filter rank=1 user_filter lxr_data.is_refseq=0

workflow annot_builder {
    take:
        gencoll_asn
        gnomon_file
        genome_asn
        parameters  // Map : extra parameter and parameter update
    main:

        def m = annot_builder_main('outdir', params).collect() 
        def i = annot_builder_input('outdir', m, '01', gnomon_file, params)
        // FIXME: intended params 4-5 to be lists of all input files and all input manifests, but it complained with only one entry
        def (all, accept, accept_ftable, annot) = annot_builder_run('outdir', i[0], gencoll_asn, i[1], gnomon_file, genome_asn, params)
    emit:
        outputs = all
        accept_asn = accept
        accept_ftable_annot = accept_ftable
        annot_files = annot
}


// [Main]
// accept_dir = "/netmnt/vast01/gpi/prod/GPIPE_PROD/data00/Lagenorhynchus_albirostris/1.1.470794/6389647/annot_builder.46361542/out/ACCEPT"
// collection_dir = "/netmnt/vast01/gpi/prod/GPIPE_PROD/data00/Lagenorhynchus_albirostris/1.1.470794/6389647/annot_builder.46361542/out/COLLECTION"
// conflict_dir = "/netmnt/vast01/gpi/prod/GPIPE_PROD/data00/Lagenorhynchus_albirostris/1.1.470794/6389647/annot_builder.46361542/out/CONFLICT"
// loss_pct_ccds = "0.0"
// loss_pct_refseq_alt_ref_loci = "100.0"
// loss_pct_refseq_other = "1.0"
// loss_pct_refseq_patches = "100.0"
// loss_pct_refseq_ref_primary = "1.0"
// reftrack_attrs_manifest = "/netmnt/vast01/gpi/prod/GPIPE_PROD/data00/Lagenorhynchus_albirostris/1.1.470794/6389647/annot_builder.46361542/inp/reftrack_attrs.mft"
// report_dir = "/netmnt/vast01/gpi/prod/GPIPE_PROD/data00/Lagenorhynchus_albirostris/1.1.470794/6389647/annot_builder.46361542/out/REPORT"
// test_dir = "/netmnt/vast01/gpi/prod/GPIPE_PROD/data00/Lagenorhynchus_albirostris/1.1.470794/6389647/annot_builder.46361542/out/TEST"

process annot_builder_main {
    input:
        val outdir
        val params
    output:
        path "annot_builder_main.ini"
    script:
    """
    #!/usr/bin/env python3
    with open('annot_builder_main.ini', 'w') as outf:
        print('[Main]', file=outf)
        print('accept_dir = "$outdir/ACCEPT"', file=outf)
        print('collection_dir = "$outdir/COLLECTION"', file=outf)
        print('conflict_dir = "$outdir/CONFLICT"', file=outf)
        print('loss_pct_ccds = "0.0"', file=outf)
        print('loss_pct_refseq_alt_ref_loci = "100.0"', file=outf)
        print('loss_pct_refseq_other = "1.0"', file=outf)
        print('loss_pct_refseq_patches = "100.0"', file=outf)
        print('loss_pct_refseq_ref_primary = "1.0"', file=outf)
        print('report_dir = "$outdir/REPORT"', file=outf)
        print('test_dir = "$outdir/TEST"', file=outf)
    """
    stub:
    """
        touch annot_builder_main.ini
        echo 'main' > annot_builder_main.ini
    """
}


// [DataProvider06]
// aliases = "Gnomon|Chainer|PartAbInitio|FullAbInitio|Chainer_GapFilled|PartAbInitio_GapFilled"
// desc = "Gnomon"
// drop_alt_brs_overlap = "1"
// enable_AR0050_AR0048 = "1"
// input_manifest = "/netmnt/vast01/gpi/prod/GPIPE_PROD/data00/Lagenorhynchus_albirostris/1.1.470794/6389647/annot_builder.46361542/inp/gnomon_models.mft"
// is_primary = "1"
// keep_top_N_models = "50"
// max_pct_ab_initio = "50"
// merge_variants = "1"
// model_maker = "gnomon2model"
// name = "gnomon"
// use_secondary_support = "1"

process annot_builder_input {
    input:
        val outdir
        path prior_file
        val provider_number
        path input_file
        val params
    output:
        path("annot_builder_input.ini") 
        path("input_manifest_${provider_number}.mft")
    script:
    """
    #!/usr/bin/env python3
    with open('annot_builder_input.ini', 'w') as outf:

        with open('${prior_file}', 'r') as f:
            print(f.read(), file=outf)

        print('[DataProvider${provider_number}]', file=outf)
        
        im = 'input_manifest_${provider_number}.mft'
        inpf = '${input_file}'
        with open(im, 'w') as mft:
            print(inpf, file=mft)
        print(f'input_manifest="{im}"', file=outf)

        print('aliases = "Gnomon|Chainer|PartAbInitio|FullAbInitio|Chainer_GapFilled|PartAbInitio_GapFilled"', file=outf)
        print('desc = "Gnomon"', file=outf)
        print('name = "gnomon"', file=outf)
        print('model_maker = "gnomon2model"', file=outf)

        print('drop_alt_brs_overlap = "1"', file=outf)
        print('enable_AR0050_AR0048 = "1"', file=outf)
        print('is_primary = "1"', file=outf)
        print('keep_top_N_models = "50"', file=outf)
        print('max_pct_ab_initio = "50"', file=outf)
        print('merge_variants = "1"', file=outf)
        print('use_secondary_support = "1"', file=outf)
    """
    stub:
    """
        touch annot_builder_input.ini
        touch input_manifest_${provider_number}.mft
        cp ${prior_file} annot_builder_input.ini
        echo 'input ${provider_number}' >> annot_builder_input.ini 
    """
}


//    ## annot_builder -accept-output both -asn-cache ${GP_cache_dir} -conffile ${conffile} -gc-assembly-manifest ${input.gencoll_asn} -logfile ${logfile}
process annot_builder_run {
    input:
        val outdir
        path conffile, stageAs: 'annot_builder_final.ini'
        path gencoll_asn
        path input_manifests
        path input_files
        path genome_asn, stageAs: 'genome/*'
        val params
    output:
        path "${outdir}/*", emit: "all"
        path "${outdir}/ACCEPT/accept.asn", emit: "accept"//, optional: true
        path "${outdir}/ACCEPT/accept.ftable_annot", emit: "accept_ftable_annot"//, optional: true
        path "${outdir}/ACCEPT/*.annot"//, optional: true
    script:
    """
    mkdir -p $outdir/ACCEPT
    mkdir -p $outdir/COLLECTION
    mkdir -p $outdir/CONFLICT
    mkdir -p $outdir/REPORT
    mkdir -p $outdir/TEST

    lds2_indexer -source genome/ -db LDS2
    # EXCEPTION_STACK_TRACE_LEVEL=Warning DEBUG_STACK_TRACE_LEVEL=Warning DIAG_POST_LEVEL=Trace
    annot_builder -accept-output both -nogenbank -lds2 LDS2 -conffile $conffile -gc-assembly $gencoll_asn -logfile ${outdir}/annot_builder.log
    cat ${outdir}/ACCEPT/*.ftable.annot > ${outdir}/ACCEPT/accept.ftable_annot
    """
    stub:
    """
    mkdir -p $outdir/ACCEPT
    mkdir -p $outdir/COLLECTION
    mkdir -p $outdir/CONFLICT
    mkdir -p $outdir/REPORT
    mkdir -p $outdir/TEST
    
    echo "1" > ${outdir}/annot_builder.log
    echo "2" > ${outdir}/accept.asn
    echo "3" > ${outdir}/accept.ftable.annot
    

    echo "4" > ${outdir}/ACCEPT/accept.asn
    echo "5" > ${outdir}/ACCEPT/accept.ftable_annot
    echo "S1" > ${outdir}/ACCEPT/S1.annot
    echo "S2" > ${outdir}/ACCEPT/S2.annot

    """
}
