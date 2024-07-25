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
