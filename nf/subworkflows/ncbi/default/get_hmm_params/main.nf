#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params } from '../../utilities'


workflow get_hmm_params {
    take:
        taxid
        parameters  // Map : extra parameter and parameter update
    main:
        def effective_taxid
        if (taxid) {
            effective_taxid = taxid
        } else {
            String params = parameters.get("gnomon_hmm_params", "")
            if (params.startsWith("taxid")) {
                effective_taxid = params.substring(6).toInteger()                
            } else {
                effective_taxid = 9606
            }
        }

        def hmm = run_get_hmm(effective_taxid)
    emit:
        outputs = hmm
}


process run_get_hmm {
    input:
        val taxid
    output:
        stdout
    script:
    """
    #!/usr/bin/env python3
    import json
    from urllib.request import urlopen
    def get_closest_hmm(taxid):
        taxon_str = str(taxid)
        if not taxon_str:
            return ""
        dataset_taxonomy_url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/"

        taxids_file = urlopen("https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/gnomon/hmm_parameters/taxid.list")
        taxids_list = []
        lineages = []
        for line in taxids_file:
            parts = line.decode("utf-8").strip().split('\t')
            if len(parts) > 0:
                t = parts[0]
                taxids_list.append(t)
                if len(parts) > 1:
                    l = map(lambda x: int(x) if x[-1] != ';' else int(x[:-1]), parts[1].split())
                    lineages.append((int(t), list(l)+[int(t)]))

        if len(lineages) < len(taxids_list):
            taxonomy_json_file = urlopen(dataset_taxonomy_url+','.join(taxids_list))
            taxonomy = json.load(taxonomy_json_file)["taxonomy_nodes"]
            lineages = [ (t["taxonomy"]["tax_id"], t["taxonomy"]["lineage"] + [t["taxonomy"]["tax_id"]]) for t in taxonomy ]

        taxon_json_file = urlopen(dataset_taxonomy_url+taxon_str)
        taxon = json.load(taxon_json_file)["taxonomy_nodes"][0]
        lineage = taxon["taxonomy"]["lineage"]
        lineage.append(taxon["taxonomy"]["tax_id"])
        # print(lineage)
        # print(taxon["taxonomy"]["organism_name"])

        best_lineage = None
        best_taxid = None
        best_score = 0
        for (t, l) in lineages:
            pos1 = 0
            last_match = 0
            for pos in range(len(lineage)):
                tax_id = lineage[pos]
                while tax_id != l[pos1]:
                    if pos1 + 1 < len(l):
                        pos1 += 1
                    else:
                        break
                if tax_id == l[pos1]:
                    last_match = pos1
                else:
                    break
            if last_match > best_score:
                best_score = last_match
                best_taxid = t
                best_lineage = l

        if best_score == 0:
            return ""
        # print(best_lineage)
        # print(best_taxid, best_score)
        return ( best_taxid,  f'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/gnomon/hmm_parameters/{best_taxid}.params')
        
    (t, url) = get_closest_hmm(${taxid})
    print(t)
    print(url)
    """
    stub:
    """
    r=\$RANDOM
    flip=\$(( r % 2 ))
    if [ "\$flip" -eq "0" ]
    then
        r=${taxid}
    fi
    touch \${r}.params
    echo \$r
    echo "file:///\${r}.params"
    """
}
