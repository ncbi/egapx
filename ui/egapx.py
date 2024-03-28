#!/usr/bin/env python
# requires pip install pyyaml

import shutil
import sys
import os
import argparse
import subprocess
import tempfile
import re
from collections import defaultdict
from pathlib import Path
from urllib.request import urlopen
import json
import yaml

VERBOSITY_DEFAULT=0
VERBOSITY_QUIET=-1
VERBOSITY_VERBOSE=1

def parse_args(argv):
    parser = argparse.ArgumentParser(description="Main script for EGAPx")
    parser.add_argument("filename", help="YAML file with input: section with at least genome: and reads: parameters")
    parser.add_argument("-e", "--executor", help="Nextflow executor, one of local, docker, aws. Uses corresponding Nextflow config file", default="local")
    parser.add_argument("-c", "--config-dir", help="Directory for executor config files, default is ./egapx_config. Can be also set as env EGAPX_CONFIG_DIR", default="")
    parser.add_argument("-o", "--output", help="Output path", default="")
    parser.add_argument("-w", "--workdir", help="Working directory for cloud executor", default="")
    parser.add_argument("-r", "--report", help="Report file prefix for report (.report.html) and timeline (.timeline.html) files, default is in output directory", default="")
    parser.add_argument("-n", "--dry-run", action="store_true", default=False)
    parser.add_argument("-q", "--quiet", dest='verbosity', action='store_const', const=VERBOSITY_QUIET, default=VERBOSITY_DEFAULT)
    parser.add_argument("-v", "--verbose", dest='verbosity', action='store_const', const=VERBOSITY_VERBOSE, default=VERBOSITY_DEFAULT)
    parser.add_argument("-fn", "--func_name", help="func_name", default="")
    return parser.parse_args(argv[1:])


def repackage_inputs(run_inputs):
    "Repackage input parameters into 'input' key if not already there"
    if 'input' in run_inputs:
        return run_inputs
    new_inputs = {}
    new_inputs['input'] = {}
    for key in run_inputs:
        if key not in { 'tasks', 'output' }:
            new_inputs['input'][key] = run_inputs[key]
        else:
            new_inputs[key] = run_inputs[key]
    return new_inputs


def convert_value(value):
    "Convert paths to absolute paths when possible, look into objects and lists"
    if isinstance(value, dict):
        return {k: convert_value(v) for k, v in value.items()}
    elif isinstance(value, list):
        return [convert_value(v) for v in value]
    else:
        if value == '' or re.match(r'^[a-z0-9]{2,5}://.*', value):
            # do not convert - this is a URL or empty string
            return value
        # convert to absolute path
        return str(Path(value).absolute())


path_inputs = { 'genome', 'hmm', 'softmask', 'reads_metadata', 'organelles', 'proteins', 'reads', 'rnaseq_alignments', 'protein_alignments' }
def convert_paths(run_inputs):
    "Convert paths to absolute paths where appropriate"
    input_root = run_inputs['input']
    for key in input_root:
        if key in path_inputs:
            if isinstance(input_root[key], list):
                input_root[key] = [convert_value(v) for v in input_root[key]]
            else:
                input_root[key] = convert_value(input_root[key])
    run_inputs['output'] = convert_value(run_inputs['output'])


def generate_reads_metadata(run_inputs):
    "Generate reads metadata file with minimal information - paired/unpaired and valid for existing libraries"
    if 'reads' not in run_inputs['input']:
        return None
    prefixes = defaultdict(list)
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
        for rf in run_inputs['input']['reads']:
            mo = re.match(r'([A-Za-z0-9]+)([^A-Za-z0-9].*)', Path(rf).parts[-1])
            if mo:
                prefixes[mo.group(1)].append(mo.group(2))
        for k, v in prefixes.items():
            paired = 'paired' if len(v) == 2 else 'unpaired'
            # SRR9005248	NA	paired	2	2	NA	NA	NA	NA	NA	NA	NA	0
            rec = "\t".join([k, 'NA', paired, '2', '2', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', '0'])
            f.write(rec + '\n')
        f.flush()
        run_inputs['input']['reads_metadata'] = f.name
        return f.name


def validate_params(run_inputs):
    "Validate input parameters"
    inputs = run_inputs['input']
    success = True
    has_genome = 'genome' in inputs
    has_rnaseq = 'reads' in inputs or 'reads_ids' in inputs or 'reads_query' in inputs
    has_hmm = 'hmm' in inputs or 'taxid' in inputs
    has_proteins = 'proteins' in inputs
    if not has_genome or not has_hmm or not (has_rnaseq or has_proteins):
        print("Missing parameter: 'genome', 'hmm', and either 'proteins', or 'reads' should be specified")
        print("  proteins and hmm can be indirectly specified by 'taxid'")
        success = False
    return success


def manage_workdir(args):
    workdir_file = f"work_dir_{args.executor}.last"
    if args.workdir:
        os.environ['NXF_WORK'] = args.workdir
        with open(workdir_file, 'w') as f:
            f.write(args.workdir)
    else:
        if os.path.exists(workdir_file):
            with open(workdir_file) as f:
                os.environ['NXF_WORK'] = f.read().strip()
        else:
            if args.executor == 'aws':
                print("Work directory not set, use -w at least once")
                return False
    return True


def get_config(script_directory, args):
    config_file = ""
    config_dir = args.config_dir if args.config_dir else os.environ.get("EGAPX_CONFIG_DIR")
    if not config_dir:
        config_dir = Path(os.getcwd()) / "egapx_config"
    if not Path(config_dir).is_dir():
        # Create directory and copy executor config files there
        from_dir = Path(script_directory) / 'assets' / 'config' / 'executor'
        if args.verbosity >= VERBOSITY_VERBOSE:
            print(f"Copy config files from {from_dir} to {config_dir}")
        if not args.dry_run:
            os.mkdir(config_dir)
            ld = os.listdir(from_dir)
            for f in ld:
                shutil.copy2(from_dir / f, config_dir)
            print(f"Edit config files in {config_dir} to reflect your actual configuration, then repeat the command")
            return ""
    config_file = Path(config_dir) / (args.executor + '.config')
    if not config_file.is_file():
        print(f"Executor {args.executor} not supported")
        return ""
    default_configs = [ "default.config" ]
    with open(str(config_file), 'r') as f:
        config_txt = f.read().replace('\n', ' ')
        # Check whether the config sets the container
        mo = re.search(r"process.+container *=", config_txt)
        if not mo:
            default_configs.append("docker_image.config")
        # Check whether the config specifies proccess memory or CPUs
        mo = re.search(r"process.+(memory|cpus) *=", config_txt)
        if not mo:
            default_configs.append("process_resources.config")
    
    # Add mandatory configs
    config_files = [str(config_file)]
    for cf in default_configs:
        config_files.append(os.path.join(script_directory, "assets/config", cf))
    return ",".join(config_files)


def get_closest_protein_bag(taxid):
    if not taxid:
        return ''
    taxon_str = str(taxid)
    dataset_taxonomy_url = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/"

    taxids_file = urlopen("https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/target_proteins/taxid.list")
    taxids_list = []
    for line in taxids_file:
        line = line.decode("utf-8").strip()
        if len(line) == 0 or line[0] == '#':
            continue
        parts = line.split('\t')
        if len(parts) > 0:
            t = parts[0]
            taxids_list.append(int(t))

    taxon_json_file = urlopen(dataset_taxonomy_url+taxon_str)
    taxon = json.load(taxon_json_file)["taxonomy_nodes"][0]
    lineage = taxon["taxonomy"]["lineage"]
    lineage.append(taxon["taxonomy"]["tax_id"])
    # print(lineage)
    # print(taxon["taxonomy"]["organism_name"])

    best_taxid = None
    best_score = 0
    for t in taxids_list:
        try:
            pos = lineage.index(t)
        except ValueError:
            continue
        if pos > best_score:
            best_score = pos
            best_taxid = t

    if best_score == 0:
        return ''
    # print(best_lineage)
    # print(best_taxid, best_score)
    return f'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/target_proteins/{best_taxid}.faa.gz'


def main(argv):
    "Main script for EGAPx"
    #warn user that this is an alpha release
    print("\n!!WARNING!!\nThis is an alpha release with limited features and organism scope to collect initial feedback on execution. Outputs are not yet complete and not intended for production use.\n")

    # Parse command line
    args = parse_args(argv)

    packaged_distro = bool(getattr(sys, '_MEIPASS', ''))
    script_directory = getattr(sys, '_MEIPASS', os.path.abspath(os.path.dirname(__file__)))
    config_file = get_config(script_directory, args)
    if not config_file:
        return 1
    
    # Check for workdir, set if not set, and manage last used workdir
    if not manage_workdir(args):
        return 1
    
    files_to_delete = []
    
    # Read default task parameters into a dict
    task_params = yaml.safe_load(open(Path(script_directory) / 'assets' / 'default_task_params.yaml', 'r'))
    run_inputs = repackage_inputs(yaml.safe_load(open(args.filename, 'r')))

    # Check for proteins input and if empty or no input at all, add closest protein bag
    if ('proteins' not in run_inputs['input'] or not run_inputs['input']['proteins']) and 'taxid' in run_inputs['input']:
        run_inputs['input']['proteins'] = get_closest_protein_bag(run_inputs['input']['taxid'])

    # Command line overrides manifest input
    if args.output:
        run_inputs['output'] = args.output
    if not validate_params(run_inputs):
        return 1

    # Convert inputs to absolute paths
    convert_paths(run_inputs)

    # Add reads_metadata.tsv file
    new_reads_metadata = generate_reads_metadata(run_inputs)
    if new_reads_metadata:
        files_to_delete.append(new_reads_metadata)

    # Add to default task parameters, if input file has some task parameters they will override the default
    task_params.update(run_inputs)

    # Move output from YAML file to arguments to have more transparent Nextflow log
    output = task_params['output']
    del task_params['output']

    if args.func_name:
        task_params['func_name'] = args.func_name

    # Run nextflow process
    if args.verbosity >= VERBOSITY_VERBOSE:
        task_params['verbose'] = True
        print("Nextflow inputs:")
        print(yaml.dump(task_params))
        if 'reads_metadata' in run_inputs['input']:
            print("Reads metadata:")
            with open(run_inputs['input']['reads_metadata'], 'r') as f:
                print(f.read())
    if packaged_distro:
        main_nf = Path(script_directory) / 'nf' / 'ui.nf'
    else:
        main_nf = Path(script_directory) / '..' / 'nf' / 'ui.nf'
    nf_cmd = ["nextflow", "-C", config_file, "-log", f"{output}/nextflow.log", "run", main_nf, "--output", output]
    if args.report:
        nf_cmd += ["-with-report", f"{args.report}.report.html", "-with-timeline", f"{args.report}.timeline.html"]
    else:
        nf_cmd += ["-with-report", f"{output}/run.report.html", "-with-timeline", f"{output}/run.timeline.html"]
    
    nf_cmd += ["-with-trace", f"{output}/run.trace.txt"]
    # if output directory does not exist, it will be created
    if not os.path.exists(output):
        os.makedirs(output)
    params_file = Path(output) / "run_params.yaml"
    nf_cmd += ["-params-file", str(params_file)]

    if args.dry_run:
        print(" ".join(map(str, nf_cmd)))
    else:
        with open(params_file, 'w') as f:
            yaml.dump(task_params, f)
            f.flush()
        if args.verbosity >= VERBOSITY_VERBOSE:
            print(" ".join(map(str, nf_cmd)))
        try:
            subprocess.run(nf_cmd, check=True, capture_output=(args.verbosity <= VERBOSITY_QUIET), text=True)
        except subprocess.CalledProcessError as e:
            print(e.stderr)
            print("To resume execution, run:")
            nf_cmd.append("-resume")
            print(" ".join(map(str, nf_cmd)))
            if files_to_delete:
                print(f"Don't forget to delete file(s) {' '.join(files_to_delete)}")
            return 1
    # TODO: Use try-finally to delete the metadata file
    for f in files_to_delete:
        os.unlink(f)

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
