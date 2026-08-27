#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { merge_params; format_multiqc_input as format_multiqc_mask_assm_stats;
    format_multiqc_input as format_multiqc_feature_counts; format_multiqc_input as format_multiqc_feature_stats;
    format_multiqc_input as format_multiqc_rnaseq_long_align_report; format_multiqc_input as format_multiqc_prot_align_stats;
    format_multiqc_input as format_multiqc_rnaseq_short_align_report; } from '../../utilities'

params.egapx_version = ""

workflow multiqc {
    take:
        star_logs
        busco_log_file
        mask_assm_stats  // Map : extra parameter and parameter update
        feature_counts_xml
        feature_stats_xml
        rnaseq_long_align_report_xml
        prot_align_stats
        rnaseq_short_align_report_xml
        parameters  // Map : extra parameter and parameter update
    main:
        String multiqc_params = merge_params('', parameters, 'multiqc')
        format_multiqc_mask_assm_stats(mask_assm_stats)
        feature_counts_xml_to_table(feature_counts_xml)
        format_multiqc_feature_counts(feature_counts_xml_to_table.out.table_xml)
        transform_feature_counts_json(format_multiqc_feature_counts.out.multiqc_inputs)
        format_multiqc_feature_stats(feature_stats_xml)
        transform_feature_stats_json(format_multiqc_feature_stats.out.multiqc_inputs)
        format_multiqc_rnaseq_long_align_report(rnaseq_long_align_report_xml)
        format_multiqc_prot_align_stats(prot_align_stats)
        format_multiqc_rnaseq_short_align_report(rnaseq_short_align_report_xml)

        get_tool_versions(params.egapx_version)
        build_software_versions_table(get_tool_versions.out.tool_versions)

        run_multiqc(star_logs.collect(), busco_log_file, format_multiqc_mask_assm_stats.out.multiqc_inputs,
                transform_feature_counts_json.out.multiqc_inputs, transform_feature_stats_json.out.multiqc_inputs,
                format_multiqc_rnaseq_long_align_report.out.multiqc_inputs, format_multiqc_prot_align_stats.out.multiqc_inputs, 
            format_multiqc_rnaseq_short_align_report.out.multiqc_inputs, build_software_versions_table.out.multiqc_inputs, multiqc_params)
    emit:
        multiqc_report = run_multiqc.out.multiqc_report
        tool_versions = get_tool_versions.out.tool_versions
}


process get_tool_versions {
    label 'single_cpu'
    label 'small_mem'
    input:
        val egapx_version
    output:
        path "output/tool_versions.yaml", emit: "tool_versions"
    script:
    """
    #!/usr/bin/env python3
    # Runs each binary's version command and parses the version string out of
    # its output, since every tool below reports its version in a different
    # format (some on stdout, some on stderr, some only via --help/-h).
    import os
    import re
    import shlex
    import subprocess

    # tool -> (command, regex with a single capture group for the version)
    commands = {
        "miniprot":    (["miniprot", "--version"], r"([\\d.]+(?:-r\\d+)?)"),
        "diamond":     (["diamond", "version"], r"diamond version\\s+([\\d.]+)"),
        "samtools":    (["samtools", "--version"], r"samtools\\s+([\\d.]+)"),
        "seqkit":      (["seqkit", "version"], r"seqkit\\s+v?([\\d.]+)"),
        "infernal":    (["cmsearch", "-h"], r"INFERNAL\\s+([\\d.]+)"),
        "tRNAscan-SE": (["tRNAscan-SE", "--version"], r"tRNAscan-SE\\s+([\\d.]+)"),
        "hmmer":       (["hmmsearch", "-h"], r"HMMER\\s+([\\d.]+)"),
        "minimap2":    (["minimap2", "--version"], r"([\\d.]+(?:-r\\d+)?)"),
        "busco":       (["busco", "--version"], r"BUSCO\\s+([\\d.]+)"),
        "STAR":        (["STAR", "--version"], r"([\\d.]+[a-z]?)"),
    }

    def run(cmd):
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            return result.stdout + result.stderr, None
        except (OSError, subprocess.TimeoutExpired) as exc:
            return "", exc

    def get_version(command, pattern):
        output, exc = run(command)
        if exc is not None:
            # Some tools (e.g. busco) are wrapper scripts that fail with
            # "Exec format error" when exec'd directly but run fine under a
            # login shell, which resolves them via PATH and sources any
            # environment setup (e.g. conda activation) they depend on.
            output, exc = run(["bash", "-lc", " ".join(shlex.quote(c) for c in command)])
        if exc is not None:
            return f"error: {exc}"

        match = re.search(pattern, output)
        return match.group(1) if match else "unknown"

    versions = {tool: get_version(command, pattern) for tool, (command, pattern) in commands.items()}
    versions["EGAPx"] = "${egapx_version}"

    os.makedirs("output", exist_ok=True)
    with open("output/tool_versions.yaml", "w") as f:
        for tool, version in versions.items():
            f.write(f"{tool}: \\"{version}\\"\\n")
    """
    stub:
    """
    mkdir -p output
    touch output/tool_versions.yaml
    """
}


process build_software_versions_table {
    label 'single_cpu'
    label 'small_mem'
    input:
        path tool_versions_yaml
    output:
        path "inputlogs/*", emit: "multiqc_inputs"
    script:
    """
    #!/usr/bin/env python3
    # Builds a MultiQC custom-content table combining each tool's version (from
    # get_tool_versions) with its Location/Authors/License/Software URL/DOI.
    # This metadata is hard-coded from build-container/Docker/LICENSE rather
    # than parsed at runtime, since that file isn't reachable from every
    # execution environment.
    #
    # Rendered as raw HTML (plot_type "html") instead of "table", since the
    # "table" plot_type always reserves a first "Sample Name" row-id column
    # that can't be removed, only relabeled.
    import json
    import os
    from html import escape

    tool_versions_file = "${tool_versions_yaml}"
    os.makedirs("inputlogs", exist_ok=True)
    outfile = os.path.join("inputlogs", "software_versions_mqc.json")

    METADATA = {
        "STAR": {
            "display_name": "STAR",
            "url": "https://github.com/alexdobin/STAR/",
            "location": "/img/gp/third-party/STAR",
            "authors": "A Dobin et al.",
            "license": "MIT License",
            "doi": "10.1093/bioinformatics/bts635",
        },
        "samtools": {
            "display_name": "SAMtools",
            "url": "https://github.com/samtools/samtools",
            "location": "/img/gp/third-party/samtools",
            "authors": "P Danecek et al.",
            "license": "The MIT/Expat License",
            "doi": "10.1093/gigascience/giab008",
        },
        "miniprot": {
            "display_name": "miniprot",
            "url": "https://github.com/lh3/miniprot",
            "location": "/img/gp/third-party/miniprot",
            "authors": "H Li",
            "license": "MIT License",
            "doi": "10.1093/bioinformatics/btad014",
        },
        "diamond": {
            "display_name": "DIAMOND",
            "url": "https://github.com/bbuchfink/diamond",
            "location": "/img/gp/third-party/diamond",
            "authors": "B Buchfink et al.",
            "license": "GPL-3.0 License",
            "doi": "10.1038/s41592-021-01101-x",
        },
        "seqkit": {
            "display_name": "SeqKit",
            "url": "https://github.com/shenwei356/seqkit",
            "location": "/img/gp/third-party/seqkit",
            "authors": "W Shen et al.",
            "license": "MIT License",
            "doi": "10.1371/journal.pone.0163962",
        },
        "infernal": {
            "display_name": "Infernal",
            "url": "http://eddylab.org/infernal/",
            "location": "/img/gp/third-party/infernal",
            "authors": "E Nawrocki et al.",
            "license": "BSD License",
            "doi": "10.1093/bioinformatics/btt509",
        },
        "tRNAscan-SE": {
            "display_name": "tRNAscan-SE",
            "url": "http://lowelab.ucsc.edu/tRNAscan-SE/",
            "location": "/img/gp/third-party/tRNAscan-SE",
            "authors": "P Chan et al.",
            "license": "GPL-3.0 License",
            "doi": "10.1093/nar/gkab688",
        },
        "hmmer": {
            "display_name": "HMMER",
            "url": "http://hmmer.org/",
            "location": "/img/gp/third-party/hmmer",
            "authors": "S Eddy",
            "license": "BSD License",
            "doi": "10.1371/journal.pcbi.1002195",
        },
        "minimap2": {
            "display_name": "minimap2",
            "url": "https://github.com/lh3/minimap2",
            "location": "/img/gp/third-party/minimap2",
            "authors": "H Li",
            "license": "MIT License",
            "doi": "10.1093/bioinformatics/bty191",
        },
        "busco": {
            "display_name": "BUSCO",
            "url": "https://busco.ezlab.org/",
            "location": "/usr/local/bin/busco",
            "authors": "Manni et al.",
            "license": "MIT License",
            "doi": "10.1093/molbev/msab199",
        },
        "EGAPx": {
            "display_name": "EGAPx",
            "url": "https://github.com/ncbi/egapx",
            "location": "https://github.com/ncbi/egapx",
            "authors": "National Center for Biotechnology Information (NCBI)",
            "license": "Public Domain (U.S. Government Work) https://github.com/ncbi/egapx/blob/main/LICENSE",
            "doi": "",
        },
    }

    def parse_tool_versions(path):
        versions = {}
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or ":" not in line:
                    continue
                tool, _, value = line.partition(":")
                versions[tool.strip()] = value.strip().strip('"')
        return versions

    def rows_to_html_table(rows, bold_columns=()):
        # Render a list of row dicts as a plain HTML <table>, with one column
        # per key (in first-seen order) and no extra row-id column.
        columns = []
        for row in rows:
            for key in row:
                if key not in columns:
                    columns.append(key)

        def render_cell(col, value):
            text = escape(str(value))
            return f"<b>{text}</b>" if col in bold_columns else text

        thead = "".join(f"<th>{escape(str(col))}</th>" for col in columns)
        tbody = "".join(
            "<tr>" + "".join(f"<td>{render_cell(col, row.get(col, ''))}</td>" for col in columns) + "</tr>"
            for row in rows
        )
        return (
            '<table class="table table-condensed table-striped mqc_table">'
            f"<thead><tr>{thead}</tr></thead>"
            f"<tbody>{tbody}</tbody>"
            "</table>"
        )

    versions = parse_tool_versions(tool_versions_file)

    rows = []
    for program, version in versions.items():
        meta = METADATA.get(program, {})
        rows.append({
            "Software Name": meta.get("display_name", program),
            "Software URL": meta.get("url", ""),
            "Location": meta.get("location", ""),
            "Authors": meta.get("authors", ""),
            "License": meta.get("license", ""),
            "Version": version,
            "DOI": meta.get("doi", ""),
        })

    out = {
        "id": "software_versions",
        "section_name": "Software Information",
        "plot_type": "html",
        "data": rows_to_html_table(rows, bold_columns={"Software Name"}),
    }
    json.dump(out, open(outfile, "w"), indent=2)
    """
    stub:
    """
    mkdir -p inputlogs
    touch inputlogs/software_versions_mqc.json
    """
}


process run_multiqc {
    label 'single_cpu'
    label 'small_mem'
    input:
        path star_logs , stageAs: 'inputlogs/*'
        path busco_log_file , stageAs: 'inputlogs/*'
        path mask_assm_stats, stageAs: 'inputlogs/*'
        path feature_counts_xml, stageAs: 'inputlogs/*'
        path feature_stats_xml, stageAs: 'inputlogs/*'
        path rnaseq_long_align_report_xml, stageAs: 'inputlogs/*'
        // Staged outside inputlogs/ so MultiQC (which only scans inputlogs/) skips it; restore stageAs: 'inputlogs/*' to re-enable this table.
        path prot_align_stats, stageAs: 'excluded_inputlogs/*'
        path rnaseq_short_align_report_xml, stageAs: 'inputlogs/*'
        path software_versions_table, stageAs: 'inputlogs/*'
        val multiqc_params
    output:
        path "output/*" , emit: 'multiqc_report'
    script:
    """
    ### Make MultiQC config file
    cat > multiqc_config.yaml << \\
    ---------------------------------------------------------------------------
    title: "EGAPx Annotation Quality Report"
    intro_text: "EGAPx is the public version of the NCBI <a href='https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/process/' target='_blank'>Eukaryotic Genome Annotation Pipeline</a>. EGAPx is available at <a href='https://github.com/ncbi/egapx' target='_blank'>GitHub</a>."

    show_analysis_paths: false
    show_analysis_time: true
    skip_generalstats: true
    skip_versions_section: true

    report_section_order:
      feature_counts_table:
        order: 8000
      feature_stats:
        order: 7000
      busco:
        order: 6000
      mask_assm_stats:
        order: 5000
      prot_align_stats:
        order: 4000
      rnaseq_align_report:
        order: 3000
      star:
        order: 2000
      long_read_report:
        order: 1000
      software_versions:
        order: -1000

    section_comments:
        feature_counts_table: "Summary counts of annotated features, which can be compared to RefSeq and other annotations of related genomes as one metric of annotation quality. Fully supported features are completely supported by experimental evidence. Features with major correction are likely protein-coding genes with corrections for premature stop codons, frameshifts, or internal gaps. A high portion of protein coding genes with major correction generally results from assembly sequence quality issues."
        feature_stats: "Summary statistics of transcript counts per gene, exon counts per transcript, and the counts and length distributions of features by subtype."
        busco: "BUSCO run in protein mode on one longest isoform per protein-coding gene as a metric of annotation quality. EGAPx automatically selects the BUSCO lineage as the odb10 dataset for the lowest taxonomic level available for the species being annotated."
        prot_align_stats: "Summary statistics of protein alignments to the genome."
        mask_assm_stats: "Percentage of the genome masked by Windowmasker."
        rnaseq_align_report: "Summary statistics of short read alignments to the genome after full EGAPx processing."
        star: "EGAPx uses custom STAR parameters which increases multimapping reads compared to default parameters. Downstream EGAPx uses custom logic to filter redundant alignments."
        long_read_report: "Summary statistics of long read alignments to the genome."

    ---------------------------------------------------------------------------

    mkdir -p output
    multiqc inputlogs $multiqc_params -o output/multiqc_report.html -n multiqc_report -o output -c multiqc_config.yaml
    """
    stub:
    """
    mkdir -p output
    mkdir -p output/multiqc_report_data
    touch output/multiqc_report.html
    touch output/multiqc_report_data/multiqc_report.txt
    """
}


process feature_counts_xml_to_table {
    label 'single_cpu'
    label 'small_mem'
    input:
        path xmlfile
    output:
        path "output/*", emit: "table_xml"
    script:
    """
    #!/usr/bin/env python3
    # Transforms a FeatureCounts XML report (as produced by EGAPx) into a flat,
    # 3-column XML table: feature_type, subtype, count
    #
    # Input XML shape:
    #   <FeatureCounts>
    #     <OverallCounts>
    #       <Counts feature_type="genes">
    #         <Total>18498</Total>
    #         <Subcount subtype="other">0</Subcount>
    #         ...
    #       </Counts>
    #       ...
    #     </OverallCounts>
    #   </FeatureCounts>
    #
    # Output XML shape:
    #   <FeatureCountsTable>
    #     <Row>
    #       <feature_type>genes</feature_type>
    #       <subtype>Total</subtype>
    #       <count>18498</count>
    #     </Row>
    #     ...
    #   </FeatureCountsTable>
    import os
    import xml.etree.ElementTree as ET

    xmlfile = "${xmlfile}"
    outname = os.path.splitext(os.path.basename(xmlfile))[0]
    os.makedirs("output", exist_ok=True)
    outfile = os.path.join("output", f"{outname}_table.xml")

    def transform(xml_path):
        # Yield (feature_type, subtype, count) rows parsed from a
        # FeatureCounts XML file.
        tree = ET.parse(xml_path)
        root = tree.getroot()

        for counts in root.iter("Counts"):
            feature_type = counts.get("feature_type", "")

            total = counts.find("Total")
            if total is not None:
                yield (feature_type, "Total", total.text)

            for subcount in counts.findall("Subcount"):
                subtype = subcount.get("subtype", "")
                count = subcount.text
                yield (feature_type, subtype, count)

    def build_table(rows):
        # Build a <FeatureCountsTable> ElementTree from
        # (feature_type, subtype, count) rows.
        root = ET.Element("FeatureCountsTable")
        for feature_type, subtype, count in rows:
            row = ET.SubElement(root, "Row")
            ET.SubElement(row, "feature_type").text = feature_type
            ET.SubElement(row, "subtype").text = subtype
            ET.SubElement(row, "count").text = count
        return ET.ElementTree(root)

    tree = build_table(transform(xmlfile))
    ET.indent(tree)
    tree.write(outfile, encoding="utf-8", xml_declaration=True)
    """
    stub:
    """
    mkdir -p output
    touch output/${xmlfile.baseName}_table.xml
    """
}


process transform_feature_counts_json {
    label 'single_cpu'
    label 'small_mem'
    input:
        path jsonfile
    output:
        path "inputlogs/*", emit: "multiqc_inputs"
    script:
    """
    #!/usr/bin/env python3
    # Transforms a MultiQC custom-content JSON file (as produced by
    # format_multiqc_input for the FeatureCounts table) so that:
    #   - the "feature_type" column is renamed to "Feature Type"
    #   - the table is rendered with ONLY the "Feature Type"/"subtype"/"count"
    #     columns -- no extra row-id column.
    #
    # NOTE: MultiQC's "table" (and "violin") plot_type always reserves and
    # renders a first column for a row/sample id -- that column can be
    # relabeled (e.g. via pconfig.col1_header) but never fully removed, since
    # every row must have a rendered identifier. The only way to drop that
    # column entirely is to stop using the "table" plot_type and instead emit
    # raw HTML via plot_type "html", where "data" is a literal HTML string
    # that MultiQC embeds as-is. That's what this process does: it builds a
    # plain <table> containing just the requested columns.
    #
    # Input "data" shape:
    #   "data": {
    #     "Row_0": {"feature_type": "genes", "subtype": "Total", "count": "18498"},
    #     ...
    #   }
    #
    # Output: the custom-content object's "plot_type" is changed to "html" and
    # "data" becomes a rendered HTML <table> string with columns "Feature
    # Type", "subtype", "count" (no row-id column). "pconfig" is dropped since
    # it only applies to the "table"/"violin" plot types.
    import json
    import os
    from html import escape

    jsonfile = "${jsonfile}"
    os.makedirs("inputlogs", exist_ok=True)
    outfile = os.path.join("inputlogs", os.path.basename(jsonfile))

    def rows_to_html_table(rows):
        # Render a list of row dicts as a plain HTML <table>, with one column
        # per key (in first-seen order) and no extra row-id column.
        columns = []
        for row in rows:
            for key in row:
                if key not in columns:
                    columns.append(key)

        thead = "".join(f"<th>{escape(str(col))}</th>" for col in columns)
        tbody = "".join(
            "<tr>" + "".join(f"<td>{escape(str(row.get(col, '')))}</td>" for col in columns) + "</tr>"
            for row in rows
        )
        return (
            '<table class="table table-condensed table-striped mqc_table">'
            f"<thead><tr>{thead}</tr></thead>"
            f"<tbody>{tbody}</tbody>"
            "</table>"
        )

    def transform(obj):
        # Return a copy of the custom-content object rendered as raw HTML
        # (plot_type "html"), with each row's "feature_type" key renamed to
        # "Feature Type" and no row-id column.
        data = obj.get("data", {})
        row_iter = data.values() if isinstance(data, dict) else data

        rows = [
            {("Feature Type" if key == "feature_type" else key): value for key, value in row.items()}
            for row in row_iter
        ]

        new_obj = dict(obj)
        new_obj.pop("pconfig", None)
        new_obj["plot_type"] = "html"
        new_obj["data"] = rows_to_html_table(rows)
        return new_obj

    with open(jsonfile) as f:
        obj = json.load(f)

    new_obj = transform(obj)

    with open(outfile, "w") as f:
        json.dump(new_obj, f, indent=2)
    """
    stub:
    """
    mkdir -p inputlogs
    touch inputlogs/${jsonfile}
    """
}


process transform_feature_stats_json {
    label 'single_cpu'
    label 'small_mem'
    input:
        path jsonfile
    output:
        path "inputlogs/*", emit: "multiqc_inputs"
    script:
    """
    #!/usr/bin/env python3
    # Transforms a MultiQC custom-content JSON file (as produced for the
    # FeatureStats table) so that:
    #   - the row-id column header ("Sample Name" by default) is renamed to
    #     "Category" via pconfig.col1_header
    #   - the "category" column is removed from every row, since it just
    #     duplicates the row-id (the dict key already used as the row label)
    #
    # Input "data" shape:
    #   "data": {
    #     "CDSs": {"category": "CDSs", "Count": "22788", "Min": "96", ...},
    #     ...
    #   }
    #
    # Output "data" shape:
    #   "data": {
    #     "CDSs": {"Count": "22788", "Min": "96", ...},
    #     ...
    #   }
    # with "pconfig.col1_header" set to "Category".
    import json
    import os

    jsonfile = "${jsonfile}"
    os.makedirs("inputlogs", exist_ok=True)
    outfile = os.path.join("inputlogs", os.path.basename(jsonfile))

    def transform(obj):
        # Return a copy of the custom-content object with the "category" key
        # stripped from every row and pconfig.col1_header set to "Category".
        data = obj.get("data", {})
        new_data = {
            row_id: {key: value for key, value in row.items() if key != "category"}
            for row_id, row in data.items()
        }

        new_obj = dict(obj)
        new_obj["data"] = new_data
        new_obj["pconfig"] = {
            **obj.get("pconfig", {}),
            "col1_header": "Category",
        }
        return new_obj

    with open(jsonfile) as f:
        obj = json.load(f)

    new_obj = transform(obj)

    with open(outfile, "w") as f:
        json.dump(new_obj, f, indent=2)
    """
    stub:
    """
    mkdir -p inputlogs
    touch inputlogs/${jsonfile}
    """
}

