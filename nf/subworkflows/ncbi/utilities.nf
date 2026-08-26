#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// Analog of shlex.split
def List<String> shellSplit(CharSequence s) {
    List<String> tokens = []
    boolean escaping = false
    char quoteChar = ' '
    boolean quoting = false
    int lastCloseQuoteIndex = Integer.MIN_VALUE
    StringBuilder current = new StringBuilder()
    
    s.eachWithIndex { c, i ->
        if (escaping) {
            current.append(c)
            escaping = false
        // } else if (c == '\\' && !(quoting && quoteChar == '\'')) {
        } else if (c == '\\' && !quoting) {
            escaping = true
        } else if (quoting && c == quoteChar) {
            quoting = false
            lastCloseQuoteIndex = i
        } else if (!quoting && (c == '\'' || c == '"')) {
            quoting = true
            quoteChar = c
        } else if (!quoting && c.isAllWhitespace()) {
            if (current.size() > 0 || lastCloseQuoteIndex == (i - 1)) {
                tokens.add(current.toString())
                current = new StringBuilder()
            }
        } else {
            current.append(c)
        }
    }
    if (current.size() > 0 || lastCloseQuoteIndex == (s.length() - 1)) {
        tokens.add(current.toString())
    }

    return tokens
}


// Convert a parameter list into a map
def Map<String, String> to_map(List<String> list )
{
    def map = [:]
    def s = list.size()
    def i = 0
    while (i < s)
    {
        def elem = list.get(i)
        i = i + 1
        if (elem.size() > 0 && elem[0] == '-')
        {
            if (i < s) {
                def val = list.get(i)
                if ( val.size() > 0 && (val[0] != '-' || val.contains(' ')) )
                {
                    map[elem] = val
                    i = i + 1
                } else {
                    map[elem] = ""
                }
            } else {
                map[elem] = ""
            }
        } else {
            println("Error: parameter string not well formed, map ${map}, elem ${elem}, i ${i}, s ${s}")
            return map
        }
    }
    return map
}


def quote(String s)
{
    if (s.size() > 0 && !(s =~ /[^\w@%+=:,.\/-]/)) {
        return s
    }
    return "'" + s + "'"
}


// Read a section of the parameters and merge them into the default parameters
// Parameters:
//   default_params: the default parameters, string
//   parameters: the parameters as a map from string to string
//   section_name: the name of the section in the parameters map to use
// Return: the merged parameters
def merge_params(default_params, parameters, section_name)
{
    def section = parameters.get(section_name, "")
    def update_map = to_map(shellSplit(section))
    def default_params_map = to_map(shellSplit(default_params))
    default_params_map.putAll(update_map)
    def l = []
    default_params_map.each { parameter, value  ->
        l << quote(parameter)
        if (value.size() > 0) {
            l << quote(value)
        }
    }

    return l.join(" ")
}



process clean_fasta_ids {
    label 'single_cpu'
    label 'small_mem'
    input:
        path fasta_in
    output:
        path "fasta_out", emit: 'fasta_out'
    script:
    """
    ## turns Fasta inputs formatted with multi-part IDs into
    ## single-part IDs, like
    ## >gi|1234|ref|NW_1234.1 Some Defline For This Org
    ## >gi|1234 Some Defline For This Org
    ## LDS chokes on the multi-part IDs.
    # the base64 nonsense is because I couldnt get it to not complain about the regex as syntax errors in some way.
    # its just this:
    # import re,sys;
    # for l in sys.stdin:
    #     <I had to delete the regex here because even in a comment nextflow lost it>
    echo 'aW1wb3J0IHJlLHN5czsKZm9yIGwgaW4gc3lzLnN0ZGluOgogICAgcHJpbnQocmUuc3ViKHIiXig+' > reol.b64
    echo 'Z2lcfFxkKylcfD8oW2Etel0rXHxbQS1aX10rW1xkXC5dK1x8KSguKikiLCAiXGc8MT5cZzwzPiIs' >> reol.b64
    echo 'IGwuc3RyaXAoKSkpCg==' >> reol.b64
    base64 -d ./reol.b64 > ./reol.py
    cat ${fasta_in} | python reol.py > ./fasta_out
    """
    stub:
    """
    touch ./fasta.out
    """
}


process multireader {
    label 'single_cpu'
    label 'small_mem'
    input:
        path fasta_file
        val parameters
    output:
        path ('output/fasta_file.asnt')  , emit: 'multireader_file'
    script:
    """
    mkdir -p output
    if [ -n "$fasta_file" ]; then
        multireader $parameters -out-format asn_text  -input $fasta_file  -output output/fasta_file.asnt
    else
         touch output/fasta_file.asnt
    fi
    """
    stub:
    """
        mkdir -p output
        touch output/fasta_file.asnt
    """
}


process convert_mask{
    label 'single_cpu'
    label 'med_mem'
    input:
        path mask
        val name
        val parameters
    output:
        path "output/${name}.asnb", emit: 'converted'
    script:
    """
    mkdir -p output
    echo "${mask}" > mask.mft
    convert_mask -input-manifest mask.mft -o output/${name}.asnb $parameters -nogenbank
    """
    stub:
    """
    mkdir -p output
    touch output/${name}.asnb
    """
}


process combine_blast_db{
    label 'single_cpu'
    label 'small_mem'
    input:
        //path contam_mask
        //path rmask_data
        //path dustmask_data
        path winmask_data
        //path rrna_mask_data
        path rfam_rrna_masks
        val name
        val parameters
    output:
        path "output/${name}.asnb", emit: 'mask_asnb'
    script:
    """
        mkdir -p output
        echo "${winmask_data}" > softmask.mft
        printf "\\n${rfam_rrna_masks}" >> softmask.mft
        combine_blast_db -input-manifest softmask.mft -o output/${name}.asnb   $parameters
    """
    stub:
    """
    mkdir -p output
    touch output/${name}.asnb
    """
}


process format_multiqc_input {
    label 'single_cpu'
    label 'small_mem'
    input:
        path xmlfile
    output:
        path "inputlogs/*", emit: "multiqc_inputs"
    script:
    """
    #!/usr/bin/env python3
    import json
    import os
    import re
    import xmltodict

    xmlfile = "${xmlfile}"
    outname = os.path.splitext(os.path.basename(xmlfile))[0]
    if not outname:
        outname = "empty_" + os.path.basename(os.getcwd())
    os.makedirs("inputlogs", exist_ok=True)
    outfile = os.path.join("inputlogs", f"{outname}_mqc.json")

    ID_ATTRS = ("@feature_type", "@subtype", "@run", "@sample", "@Taxid", "@taxid",
                "@category", "@Category", "@Name", "@name",
                "@accession", "@Accession")

    def strip_attr_prefix(key):
        return key[1:] if key.startswith("@") else key

    def flatten(obj, prefix=""):
        # Recursively flatten a nested dict/list produced by xmltodict into a
        # single-level dict of "dotted.key" -> scalar value, for use as the
        # columns of a single table row. Nested lists here are folded into
        # columns (not new rows) -- row-splitting is handled by extract_rows().
        flat = {}
        if isinstance(obj, dict):
            for key, value in obj.items():
                if key == "#text":
                    flat[prefix or "value"] = value
                    continue
                child_key = strip_attr_prefix(key)
                new_prefix = f"{prefix}.{child_key}" if prefix else child_key
                flat.update(flatten(value, new_prefix))
        elif isinstance(obj, list):
            for i, item in enumerate(obj):
                label = None
                if isinstance(item, dict):
                    for id_attr in ID_ATTRS:
                        if id_attr in item:
                            label = item[id_attr]
                            break
                if label is None:
                    label = str(i)
                new_prefix = f"{prefix}.{label}" if prefix else label
                flat.update(flatten(item, new_prefix))
        else:
            flat[prefix or "value"] = obj
        return flat

    def contains_list(node):
        # True if a list appears anywhere within node (distinguishes a
        # "container" element from a "leaf" record that can be flattened
        # as-is into a single row).
        if isinstance(node, list):
            return True
        if isinstance(node, dict):
            return any(contains_list(v) for k, v in node.items() if k != "#text")
        return False

    def pick_row_label(item, fallback):
        if isinstance(item, dict):
            for id_attr in ID_ATTRS:
                if id_attr in item:
                    return str(item[id_attr])
        return str(fallback)

    def unique_row_label(label, rows):
        row_label = label
        suffix = 2
        while row_label in rows:
            row_label = f"{label} ({suffix})"
            suffix += 1
        return row_label

    def extract_rows(node, container_name):
        # Recursively walk down through "container" elements until repeated
        # (list) XML elements are found; each item of such a list becomes
        # its own row in the output table, so that multiple rows in the XML
        # become multiple rows in JSON, instead of being merged into one.
        rows = {}
        loose = {}
        if not isinstance(node, dict):
            return rows

        for key, value in node.items():
            if key == "#text":
                loose[container_name] = value
                continue
            col_key = strip_attr_prefix(key)
            if isinstance(value, list):
                for i, item in enumerate(value):
                    label = unique_row_label(pick_row_label(item, f"{col_key}_{i}"), rows)
                    rows[label] = flatten(item)
            elif isinstance(value, dict):
                if contains_list(value):
                    rows.update(extract_rows(value, col_key))
                else:
                    rows[col_key] = flatten(value)
            else:
                loose[col_key] = value

        if loose:
            if len(rows) == 1:
                next(iter(rows.values())).update(loose)
            else:
                rows[unique_row_label(container_name, rows)] = loose

        return rows

    try:
        parsed_dict = xmltodict.parse(open(xmlfile).read())
    except Exception as e:
        print(f"Error parsing XML file {xmlfile}: {e}")
        json.dump({}, open(outfile, "w"))
        exit(0)

    if not parsed_dict:
        json.dump({}, open(outfile, "w"))
        exit(0)

    root_key = next(iter(parsed_dict))
    root_content = parsed_dict[root_key]

    def camel_to_title(name):
        # "FeatureCountsTable" -> "Feature Counts Table"
        return re.sub(r'(?<!^)(?=[A-Z])', ' ', name)

    rows = extract_rows(root_content, root_key)
    if not rows:
        rows = {outname: flatten(root_content)}
    elif len(rows) == 1:
        rows = {outname: next(iter(rows.values()))}

    out = {
        "id": outname,
        "section_name": camel_to_title(root_key),
        "plot_type": "table",
        "pconfig": {
            "id": f"{outname}_table",
            "namespace": root_key,
        },
        "data": rows,
    }
    json.dump(out, open(outfile, "w"), indent=2)
    """
}
