#!/usr/bin/env nextflow
nextflow.enable.dsl=2


include {   checkpoint_save; 
            checkpoint_load_channels; 
            checkpoint_load_structure;
            } from './checkpoint_both.nf'

include { dummy_checkpoint_source; 
         complicated_dummy_checkpoint_source;
         complicated_dummy_checkpoint_source_cached;
         dummy_checkpoint_source_cached;
         ECHO;   
         ECHO as ECHO_01; ECHO as ECHO_02;    
          } from './checkpoint_test_wfs.nf'


workflow test_checkpoint_save {
    main:
    def ckName = 'save_test'
    def ckDir = "${projectDir}/.save_test_ck"
    
    dummy_checkpoint_source()
    checkpoint_save([dummy_checkpoint_source.out], [name: ckName, dir: ckDir])
    
    // Verify files were written
    checkpoint_save.out.saved_mft.view { "Saved file: $it" }
    
    emit:
    saved_mft = checkpoint_save.out.saved_mft
}

workflow test_checkpoint_save_complicated {
 main:
    def ckName = 'save_test'
    def ckDir = "${projectDir}/.save_test_ck"
    
    complicated_dummy_checkpoint_source()
    checkpoint_save([complicated_dummy_checkpoint_source.out], [name: ckName, dir: ckDir])
    
    // Verify files were written
    checkpoint_save.out.saved_mft.view { "Saved file: $it" }
    
    emit:
    saved_mft = checkpoint_save.out.saved_mft   
}

workflow test_checkpoint_load {
    main:
    def ckName = 'save_test'  // Must match the save test above
    def ckDir = "${projectDir}/.save_test_ck"
    
    loaded_struct = checkpoint_load_structure([name: ckName, dir: ckDir])
    //println(loaded_struct)
    assert loaded_struct && !loaded_struct.isEmpty() : "Load returned empty map"
    assert loaded_struct.path_single_q.kind == 'queue'
    assert loaded_struct.path_list_q.kind   == 'queue'
    assert loaded_struct.path_single_v.kind == 'value'
    assert loaded_struct.path_list_v.kind   == 'value'

    loaded_channels = checkpoint_load_channels([name: ckName, dir: ckDir])
    per_key = loaded_channels.collect { k, ch ->
        ch.ifEmpty([]).toList().map { vals ->
            [[k, vals.size()]]
        }
    }
    merged = per_key.inject(channel.empty()) { acc, c -> acc.mix(c) }
    summary = merged.collect().map { rows ->
        //println(rows)
        def seen = rows.collect { it[0] } as Set
        assert seen as Set == loaded_channels.keySet() as Set :
            "Loaded key mismatch seen=${seen.sort()} expected=${loaded_channels.keySet().sort()}"
        def dupKeys = rows.groupBy { it[0] }.findAll { k,v -> v.size() > 1 }.keySet()
        assert dupKeys.isEmpty() : "Duplicate key emissions: ${dupKeys}"
        rows
    }
    summary.view { rows -> "Load validation ok: ${rows}" }

    emit:
    out = loaded_channels
 }

// rm the cache dir, 
// run this twice
// should just be identical results, and file exists
workflow test_checkpoint_cached {
    main:
    def ckDir = "${projectDir}/.cached_test_ck"
    
    dummy_checkpoint_source_cached(ckDir)

    dummy_checkpoint_source_cached.out.strings.view { "CACHED output: $it" }

    ECHO(dummy_checkpoint_source_cached.out.strings)


    dummy_checkpoint_source_cached.out.null_six.view { "CACHED six output: $it" }


    emit:
    out = ECHO.out.echoed
}

workflow test_checkpoint_cached_complicated {
    main:
    def ckDir = "${projectDir}/.cached_test_ck"
    
    complicated_dummy_checkpoint_source_cached(ckDir)

    complicated_dummy_checkpoint_source_cached.out.letters.view { "CACHED output: $it" }

    ECHO(complicated_dummy_checkpoint_source_cached.out.letters)

    emit:
    out = ECHO.out.echoed
}

process RM {
    input:
    val dir

    output:
    val dir , emit: dir_out 

    script:
    """
    rm -I -rf ${dir}
    mkdir -p ${dir}
    """
    
}

workflow passthrough {
    take: ch
    main:
    out = ch.collect().flatMap { it }
    emit:
    out = out
}


// this is the closest it gets to a 'round-trip test'
workflow test_save_then_load {
    main:
    def ckName = 'roundtrip_test'
    def ckDir  = "${projectDir}/.roundtrip_ck"

    RM(ckDir)
    dummy_checkpoint_source()
    checkpoint_save([dummy_checkpoint_source.out], [name: ckName, dir: ckDir])

    // gate everything on save completion
    def ready = checkpoint_save.out.saved_mft.collect()

    // load structure (plain function, allowed in closure)
    ld = checkpoint_load_channels([name: ckName, dir: ckDir])
    def loaded_ch = ld

    println("Workflow structure: ${dummy_checkpoint_source.out.channels.keySet()}")
    loaded_ch.view { "Loaded structure: ${it.keySet()}" }
    // test one channels
    def path_list_ch    = loaded_ch.flatMap { s -> s['path_list'] }

    ECHO_01(dummy_checkpoint_source.out.path_list)
    ECHO_02(path_list_ch)
    
    def filenames = { ch ->
        ch
        .flatMap { v -> (v instanceof Collection) ? v : [v] }   // normalize list-vs-item
        .map { p -> file(p).getName() }                          // true basename
        .collect()
        .map { [it.sort()] }
        }

    src_names_ch    = filenames(dummy_checkpoint_source.out.path_list)
    loaded_names_ch = filenames(path_list_ch)
    src_names_ch
        .combine(loaded_names_ch)
        .map { left, right ->
            def src = new ArrayList(left).sort()
            def loaded = new ArrayList(right).sort()
            assert src == loaded : "Mismatch src=${src} loaded=${loaded}"
            "path_list ok: ${src}"
        }
        .view()

    emit:
    saved = checkpoint_save.out.saved_mft
}

workflow test_all {
    main:
    test_checkpoint_save()
    test_checkpoint_save_complicated()

    test_checkpoint_load()
    test_checkpoint_cached()
    // cant just call it twice
    //test_checkpoint_cached()
    test_save_then_load()
}
///