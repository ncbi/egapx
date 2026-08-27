nextflow.enable.dsl=2
import groovy.json.JsonOutput
import groovy.json.JsonSlurper


import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.DataflowQueue
import groovyx.gpars.dataflow.DataflowBroadcast


def to_json_safe(obj) {
    switch (obj) {
        case null:
            return null
        case Number:
        case Boolean:
        case CharSequence:
            return obj
        case java.nio.file.Path:
            return [ __type__: 'path', value: obj.toString() ]
        case File:
            return [ __type__: 'path', value: obj.toString() ]
        case Tuple:
            return [ __type__: 'tuple', value: obj.toList().collect { to_json_safe(it) } ]
        case Map:
            return obj.collectEntries { k, v -> [ k.toString(), to_json_safe(v) ] }
        case Iterable:
            return obj.collect { to_json_safe(it) }
        default:
            return obj.toString()
    }
}

def from_json_safe(obj) {
    switch (obj) {
        case null:
            return null
        case Number:
        case Boolean:
        case CharSequence:
            return obj
        case Map:
            if (obj.__type__ == 'path')
                return java.nio.file.Paths.get(obj.value.toString())
            if (obj.__type__ == 'tuple')
                return tuple(*((obj.value as List).collect { from_json_safe(it) }))
            return obj.collectEntries { k, v -> [ k, from_json_safe(v) ] }
        case List:
            return obj.collect { from_json_safe(it) }
        default:
            return obj
    }
}


boolean is_dataflow(obj) {
    obj instanceof DataflowReadChannel ||
    obj instanceof DataflowVariable   ||
    obj instanceof DataflowQueue      ||
    obj instanceof DataflowBroadcast
}

def infer_channel_kind(obj) {
    //println("ick: ${obj} : ${obj.getClass()}")
    if (obj == null) return 'none'                 
    if (obj instanceof DataflowVariable) {
        //println("is dataflow variable, bound=${obj.isBound()}, length ${obj.metaClass}")
        return 'value'
    }
    if (is_dataflow(obj)) return 'queue'
    return 'value'                                  
}

def to_channel(obj) {
    //println("${obj}, ${obj.getClass()} , ${obj.dump()}")

    if (obj instanceof Channel)      return obj          // already a channel
    if (obj == null)                 return channel.empty()
    if (is_dataflow(obj))             return obj   
    if (obj instanceof Collection ||
        obj.getClass().isArray())    return channel.value(obj)  // expand collection
    return channel.of(obj)                                   // treat as scalar
}

/* ----------------------  SAVE PROCESS  ---------------------- */

process CHECKPOINT_SAVE {
    publishDir "${dir}", pattern: "*.json"
    
    input:
    tuple val(name), val(json_text)
    val name_arg
    val dir

    output:
    path("${name_arg}.json"), emit: saved_file

    script: 
    """
    printf '%s\\n' '${json_text}' > ${name_arg}.json
    """
    
    //exec:
    //    file("${name_arg}.json").toFile().setText(json_text + "\n")

}

/* ----------------------  WORKFLOW: checkpoint_save  ---------------------- */

workflow checkpoint_save {
    take: wf_out
    take: chpt_args

    main:
    def returned_file = null
    if (!chpt_args.enabled) {
        returned_file = null   
    } else {
        Map norm = Checkpoint.normalizeOut(wf_out)
        def orig_keys = norm.keySet() as List

        def static_rows = []
        def channel_rows = []
        def channel_kinds = [:]

        // Separate already resolved values from dynamic channels
        // for things like real constants, that do not have to flow through anything
        norm.each { key, v ->
            def kind = infer_channel_kind(v)
            channel_kinds[key] = kind

            // Preserve "plain null" as non-channel
            if (kind == 'none') {
                return
            }
            // Bound value channels are safe to read directly
            if (v instanceof DataflowVariable && v.isBound()) {
                
                //    println("${v}, ${v.getClass()}")
                //    println("decl: ${v.declaredClass.declaredMethods}")
                
                static_rows << [(key): [__item__: to_json_safe(v.getVal())]]
                return
            }
            // Everything else flows through channel plumbing
            def ch = to_channel(v)
            channel_rows << ch.map { item ->
                //println("ctoc: ${item} : ${item.getClass()}")
                [(key): [__item__: to_json_safe(item)]]
            }
        }

        def static_ch = static_rows ? channel.from(static_rows) : channel.empty()
        def stream_ch = channel_rows.inject(channel.empty()) { acc, ch -> acc.mix(ch) }

        def merged_ch = static_ch.mix(stream_ch)

        // channel.empty() values will have disappeared from here, 
        // so the output map is seeded with the original input keys, so anything that
        // isnt repopulated by the channels will be preserved in the output. 

        def json_ch = merged_ch.collect().map { rows ->
            //def json_root = orig_keys.collectEntries { k -> [(k): []] }
            def observed = orig_keys.collectEntries { k ->
                [(k): [kind: channel_kinds[k], seen: false, items: []]]
            }

            rows.each { row ->
                //late_kind = infer_channel_kind(row)
                //println("late kind: ${late_kind} for rows: ${row}")
                
                row.each { k, payload ->
                    observed[k].seen = true
                    observed[k].items << payload.__item__
                }
            }
            
            // Final serialized shape per output key
            def channels_obj = [:]
            orig_keys.each { k ->
                def e = observed[k]
                def kind = e.kind
                def seen = e.seen
                def items = e.items

                switch (kind) {
                    case 'none':
                        channels_obj[k] = [kind: 'none', seen: seen]
                        break

                    case 'value':
                        // Preserve value semantics, including null value
                        def value = seen ? (items.size() == 1 ? items[0] : items) : null
                        channels_obj[k] = [kind: 'value', seen: seen, value: value]
                        break

                    case 'queue':
                        // No emissions observed => true empty channel
                        if (!seen) {
                            channels_obj[k] = [kind: 'empty', seen: seen]
                        } else {
                            channels_obj[k] = [kind: 'queue', seen: seen, items: items]
                        }
                        break

                    default:
                        throw new IllegalStateException("Unsupported channel kind '${kind}' for key '${k}'")
                }
            }
            def json_text = JsonOutput.prettyPrint(JsonOutput.toJson([
                name: chpt_args.name,
                channels: channels_obj,
                dir : chpt_args.dir
            ]))
            tuple(chpt_args.name, json_text)
        }

        CHECKPOINT_SAVE(json_ch, chpt_args.name, chpt_args.dir)
        returned_file = CHECKPOINT_SAVE.out.saved_file
    }

    emit:
    saved_mft = returned_file
}

/* ----------------------  WORKFLOW: checkpoint_load  ---------------------- */
/*
def reconstruct_channel(payload) {
    if (payload == null) return channel.empty()

    if (payload instanceof Map && payload.containsKey('kind') && payload.containsKey('items')) {
        def kind = payload.kind?.toString() ?: 'queue'
        def items = (payload.items instanceof List) ? payload.items : [payload.items]
        def restored = items.collect { from_json_safe(it) }
        switch (kind) {
            case 'empty':
                return channel.empty()
            case 'value':
                return restored.size() == 1 ? channel.value(restored[0]) : channel.value(restored)
            default:
                return restored.size() == 1 ? channel.of(restored[0]) : channel.from(restored)
        }
    }

    def list = (payload instanceof List) ? payload : [payload]
    def restored = list.collect { from_json_safe(it) }
    return restored.size() == 1 ? channel.value(restored[0]) : channel.from(restored)
}
*/
def reconstruct_channel(payload) {
    switch (payload.kind) {
        case "empty":
            return channel.empty()
        case "queue":
            return payload.items.size() == 1 ? channel.of(payload.items[0]) : channel.from(payload.items)
        case "value":
            if(!payload.seen) return channel.empty()
            return channel.value(payload.value)
        case "scalar":
            return payload.value // not a channel, just a value
        case "none":
            return null
        default:
            throw new IllegalArgumentException("Unknown channel kind: ${payload.kind}")
    }
}

def infer_legacy_kind(vals) {
    def list = (vals instanceof List) ? vals : [vals]
    return (list.size() <= 1) ? 'value' : 'queue'
}


// Channel-facing API
def checkpoint_load_channels(chpt_args) {
    def structure = checkpoint_load_structure(chpt_args)
    return structure.collectEntries { k, v ->
        [(k): reconstruct_channel(v)]
    }
}

def checkpoint_load_structure(chpt_args) {
    def pattern = "${chpt_args.dir}/${chpt_args.name}.json"
   
    if (file(pattern).exists()) {
        def parsed = new JsonSlurper().parseText(file(pattern).text)
        def root = parsed.channels 
        if (!(root instanceof Map)) {
            throw new IllegalStateException("Invalid checkpoint: missing channels map in ${pattern}")
        }
        return from_json_safe(root)
    }

    return [:]
}

process CHECKPOINT_LOAD_JSON {
    tag "${name_arg}"
    input:
    val name_arg
    val dir_arg

    output:
    path("${name_arg}.json"), emit: json_file

    script:
    """
    set -euo pipefail
    src="${dir_arg}/${name_arg}.json"
    [[ -f "\$src" ]] || { echo "Missing checkpoint: \$src" >&2; exit 1; }
    cp "\$src" "${name_arg}.json"
    """
}


///
/*
this does not work
make a 
workflow_name_cached workflow
workflow checkpoint_load {
    take:
    chpt_args

    main:
    CHECKPOINT_LOAD_JSON(chpt_args.name, chpt_args.dir)

    // hide channel plumbing here
    loaded = CHECKPOINT_LOAD_JSON.out.json_file
        .map { jf -> checkpoint_load_channels(jf) }
        .first()

    emit:
    loaded = loaded
}
*/
/* ---------------------- END FILE ---------------------- */
