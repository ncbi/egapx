#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {   checkpoint_save; 
            checkpoint_load_channels; 
            checkpoint_load_structure;
            CHECKPOINT_LOAD_JSON} from './checkpoint_both.nf'


process ECHO {
    input:
        val(msg)
    output:        
        stdout emit: echoed
    script:
    """
    echo ECHO "${msg}"
    """
}

process TEST_CHECKPOINT_SAVE {
    input:
    output:
        path('single.yaml'), emit: 'path_single'
        path("*.txt"), emit: 'path_list'    
    script:
    """
    echo 'a' > a.txt
    echo 'b' > b.txt
    echo 'c' > c.txt
    echo 'd' > d.txt
    echo 'single.yaml' > single.yaml
    """
}

process TEST_CHECKPOINT_MANY_A {
    input:
    output:
        path("*.txt"), emit: 'path_list'    
    script:
    """
    echo 'a' > a.txt
    echo 'b' > b.txt
    echo 'c' > c.txt
    echo 'd' > d.txt
    echo 'e' > e.txt
    echo 'f' > f.txt
    echo 'g' > g.txt
    echo 'h' > h.txt
    """
}

process TEST_CHECKPOINT_MANY_B {
    input:
        val letter
    output:
        path("*.txt"), emit: 'path_list'    
    script:
    """
    echo '${letter}' > ${letter}.txt
    """
}

process TEST_CHECKPOINT_ONE {
    input:
    output:
        path('single.yaml'), emit: 'path_single'
    script:
    """
    echo 'single.yaml' > single.yaml
    """
}

process TEST_CHECKPOINT_USE_BOTH {
    input:
    path letter_file
    path yaml_file
    output:
      path("out/*.txt"), emit: 'letter_yaml'
    script:
    """
    mkdir -p out
    echo "letter file: ${letter_file}, yaml file: ${yaml_file}" > out/result.txt
    """
}

process REPEATER {
    input:
        path msg_file
    output:
        val msg_file, emit: repeated
    script:
    """
    touch ${msg_file}
    """
}

workflow dummy_checkpoint_source {
    main:
    TEST_CHECKPOINT_SAVE()
    
    def ch_path_single_q = TEST_CHECKPOINT_SAVE.out.path_single
    def ch_path_list_q   = TEST_CHECKPOINT_SAVE.out.path_list

    // explicit value semantics derived from queue channels
    def ch_path_single_v = ch_path_single_q.first()   // one Path value
    def ch_path_list_v   = ch_path_list_q.toList()    // List<Path> value

    def y = EMPTY_LIST()
    def x = never_runs(y.out.empty_file_list)

    emit:
    empty = channel.empty()
    null_item = channel.of(null)
    ///null_two  = channel.value(null)       // hang
    null_thre = channel.of([null,null,null])
    null_four = channel.value([null,null,null])  
    null_five = channel.of(null, null, null)
    null_six  = null
    null_sevn = channel.value([null])
    numbers   = channel.of(10, 20, 30)
    strings   = ['apple', 'orange', 'banana']
    number    = 42
    diff_nums = [1, 2.001, 'three']
    mixed     = channel.of(1, 2.001, 'three')
    deep      = channel.of([a:[b:[c:1]], list:[[1,2],[3,4]]])
    dup_keys  = channel.of([k:1], [k:2])
    tuple_w_path = channel.of(tuple('id1', [meta:[x:1], vals:[1,2], p:file('a.txt')]))
    text      = channel.of('hello')
    tuples    = channel.of( tuple('a',1), tuple('b',2) )
    mapval    = channel.of( [x:1, y:2] )
     // keep old names if you want existing tests unchanged
    //path_single = path_single_q
    //path_list   = path_list_q
    // new explicit compare pair
    path_single_q = ch_path_single_q
    path_single_v = ch_path_single_v
    path_list_q   = ch_path_list_q
    path_list_v   = ch_path_list_v
    never_runs    = x
}

process run_align_report {
    input:
        path gencoll_asn
        path input_metadata
        path intron_counts
        path run_stats
        val run_list
        val params
    output:
        path "rnaseq_align_report.xml", emit: "align_report"
        path "*_runs.txt", emit: "run_reports"
    script:
    """
    touch rnaseq_align_report.xml
    touch stub_report_runs.txt
    """
}

process EMPTY_LIST {
    output:
    path '*.txt', emit: empty_file_list, arity:0..100
    script:
    """
    rm *.txt  || true
    """
}

workflow never_runs {
    take:
        x
    main:
        z = x.count
    emit:
        outfiles = z
}

workflow dummy_checkpoint_source_cached {
    take:
    checkpoint_dir

    main:
    def save_name = 'chpt_test'
    def restored

    if (file("${checkpoint_dir}/${save_name}.json").exists()) {
        //checkpoint_load([name: save_name, dir: checkpoint_dir])
        //def loaded_ch = checkpoint_load.out
        //restored = loaded_ch
        CHECKPOINT_LOAD_JSON(save_name, checkpoint_dir)
        println("cache hit: ${checkpoint_dir}")
        restored = checkpoint_load_channels([name: save_name, dir: checkpoint_dir])
    }
    else {
        println("cache miss: ${checkpoint_dir} -> running source + save")
        dummy_checkpoint_source()
        checkpoint_save([dummy_checkpoint_source.out], [name: save_name, dir: checkpoint_dir])
        restored = Checkpoint.normalizeOut(dummy_checkpoint_source.out)
    }


    emit:
    empty        = restored.empty
    null_item    = restored.null_item
    ///null_two     = restored.null_two
    null_thre    = restored.null_thre 
    null_four    = restored.null_four
    null_five    = restored.null_five
    null_six     = restored.null_six
    null_sevn    = restored.null_sevn
    numbers      = restored.numbers
    strings      = restored.strings
    number       = restored.number
    diff_nums    = restored.diff_nums
    mixed        = restored.mixed
    deep         = restored.deep
    dup_keys     = restored.dup_keys
    tuple_w_path = restored.tuple_w_path
    text         = restored.text
    tuples       = restored.tuples
    mapval       = restored.mapval
    path_single  = restored.path_single
    path_list    = restored.path_list
    never_runs   = restored.never_runs
}


workflow complicated_dummy_checkpoint_source {
    main:
    /////TEST_CHECKPOINT_SAVE()
    
    TEST_CHECKPOINT_MANY_A()
    def ma = TEST_CHECKPOINT_MANY_A.out.path_list
    //TEST_CHECKPOINT_MANY_B(channel.of('i','j','k','l','m','n','o','p'))
    //def ma = TEST_CHECKPOINT_MANY_B.out.path_list
    TEST_CHECKPOINT_ONE()
    //TEST_CHECKPOINT_USE_BOTH(TEST_CHECKPOINT_MANY.out.path_list.flatten(), TEST_CHECKPOINT_ONE.out.path_single)

    //println(TEST_CHECKPOINT_MANY.out.many)
    //REPEATER(TEST_CHECKPOINT_MANY.out.many)
    
    emit:
    many = ma
    one  = TEST_CHECKPOINT_ONE.out.path_single
    //rr   = REPEATER.out.repeated
    //complicated_path_list = TEST_CHECKPOINT_USE_BOTH.out.letter_yaml
   

}

workflow complicated_dummy_checkpoint_source_cached {
    take:
    checkpoint_dir

    main:
    def save_name = 'comp_chpt_test'
    def restored

    if (file("${checkpoint_dir}/${save_name}.json").exists()) {
        CHECKPOINT_LOAD_JSON(save_name, checkpoint_dir)
        println("cache hit: ${checkpoint_dir}")
        restored = checkpoint_load_channels([name: save_name, dir: checkpoint_dir])
    }
    else {
        println("cache miss: ${checkpoint_dir} -> running source + save")
        complicated_dummy_checkpoint_source()
        checkpoint_save([complicated_dummy_checkpoint_source.out], [name: save_name, dir: checkpoint_dir])
        restored = Checkpoint.normalizeOut(complicated_dummy_checkpoint_source.out)
    }
    
    //def rr = REPEATER(restored.many))

    TEST_CHECKPOINT_USE_BOTH(restored.many, restored.one)
    //TEST_CHECKPOINT_USE_BOTH(restored.many.flatten(), restored.one)

    //checkpoint_save([TEST_CHECKPOINT_USE_BOTH.out], [name: 'comp_both_test', dir: checkpoint_dir])

    emit:    
        many  = restored.many
        one   = restored.one
        letters = TEST_CHECKPOINT_USE_BOTH.out.letter_yaml
}
//