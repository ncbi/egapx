#!/usr/bin/env python
# requires pip install pyyaml.
#
from errno import ENOENT
import shlex
import shutil
import sys
import os
import argparse
import subprocess
import re
import time
import datetime
from collections import defaultdict
from ftplib import FTP
import ftplib
from pathlib import Path
from typing import List
from urllib.request import urlopen
from urllib.parse import quote
import json
import csv
import sqlite3
import glob
import hashlib
import tarfile
import urllib.request

# Requires pip install pyyaml
import yaml

software_version = "0.4.0-alpha"

start_time = time.time()

VERBOSITY_DEFAULT=0
VERBOSITY_QUIET=-1
VERBOSITY_VERBOSE=1

# Clades with special treatment
MAGNOLIOPSIDA = 3398
ACTINOPTERYGII = 7898
COELACANTHIMORPHA = 118072
CHONDRICHTHYES = 7777
DIPNOMORPHA = 7878
FISH = {ACTINOPTERYGII, COELACANTHIMORPHA, CHONDRICHTHYES, DIPNOMORPHA}
PRIMATE = 9443
MAMMAL = 40674
RODENT = 9989
VERTEBRATE = 7742
MAGNOLIOPSIDA = 3398
ANIMALS = 33208   
## define gbdiv_inv as has ANIMALS but not VERTEBRATE
INSECTA = 50557
ARTHROPOD = 6656
VIRIDIPLANTAE = 33090
LEPIDOSAURS = 8504 
AMPHIBIANS = 8292 
ECHINODERMATA = 7586

BUSCO_DATA_URL = "https://busco-data.ezlab.org/v5/data"
FTP_EGAP_PROTOCOL = "https"
FTP_EGAP_SERVER = "ftp.ncbi.nlm.nih.gov"
FTP_EGAP_ROOT_PATH = "genomes/TOOLS/EGAP/support_data"
FTP_EGAP_ROOT = f"{FTP_EGAP_PROTOCOL}://{FTP_EGAP_SERVER}/{FTP_EGAP_ROOT_PATH}"
DATA_VERSION = "current"
dataset_taxonomy_url = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/"
SRA_DOWNLOAD_FOLDER = 'sra_dir'
SRA_RUNS_FILE = 'runs.yaml'
SRA_METADATA_FILE = 'metadata.yaml'
SRA_QUERIES_FILE = 'queries.yaml'

# Soft limit for SRA queries, can be overridden by --force
RNASEQ_QUERY_LIMIT = int(os.environ.get('EGAPX_RNASEQ_QUERY_LIMIT', '20'))

esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&retmode=json&retmax=2147483647&idtype=gi&term="
runinfo_url = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?retmode=csv&sp=runinfo&uid="
TAB=chr(9)
NL=chr(10)

user_cache_dir = ''


def parse_args(argv):
    parser = argparse.ArgumentParser(description="Main script for EGAPx")
    group = parser.add_argument_group('run')
    group.add_argument("filename", nargs='?', help="YAML file with input: section")
    group.add_argument("-o", "--output", help="Output path", default="")
    parser.add_argument("-e", "--executor", help="Nextflow executor, one of docker, singularity, aws, or local (for NCBI internal use only). Uses corresponding Nextflow config file", default="local")
    parser.add_argument("-c", "--config-dir", help="Directory for executor config files, default is ./egapx_config. Can be also set as env EGAPX_CONFIG_DIR", default="")
    parser.add_argument("-w", "--workdir", help="Working directory for cloud executor", default="")
    parser.add_argument("-r", "--report", help="Report file prefix for report (.report.html) and timeline (.timeline.html) files, default is in output directory", default="")
    parser.add_argument("-n", "--dry-run", action="store_true", default=False)
    parser.add_argument("-st", "--stub-run", action="store_true", default=False)
    parser.add_argument("-so", "--summary-only", help="Print result statistics only if available, do not compute result", action="store_true", default=False)
    parser.add_argument("-lo", "--logs-only", help="Collect execution logs if available, put them in output directory, do not compute result", action="store_true", default=False)
    parser.add_argument("-ot", "--ortho-taxid", default=0, help="Taxid of reference data for orthology tasks")
    group = parser.add_argument_group('download')
    group.add_argument("-dl", "--download-only", help="Download external files to local storage, so that future runs can be isolated", action="store_true", default=False)
    parser.add_argument("-lc", "--local-cache", help="Where to store the downloaded files", default="")
    parser.add_argument("-q", "--quiet", dest='verbosity', action='store_const', const=VERBOSITY_QUIET, default=VERBOSITY_DEFAULT)
    parser.add_argument("-v", "--verbose", dest='verbosity', action='store_const', const=VERBOSITY_VERBOSE, default=VERBOSITY_DEFAULT)
    parser.add_argument("-V", "--version", dest='report_version', help="Report software version", action='store_true', default=False)
    parser.add_argument("-fn", "--func_name", help="func_name", default="")
    parser.add_argument("--force", dest='force', help="Force execution despite of warnings", action='store_true', default=False)
    
    return parser.parse_args(argv[1:])


class FtpDownloader:
    def __init__(self):
        self.host = None
        self.ftp = None

    def connect(self, host):
        self.host = host
        self.reconnect() 

    def reconnect(self):
        self.ftp = FTP(self.host)
        self.ftp.login()
        self.ftp.set_debuglevel(0)
       
    ##ftp_types = set()
    def download_ftp_file(self, ftp_name, local_path):
        # print(f"file: {ftp_name}")
        # print(f"f: { os.path.dirname(local_path)}")

        os.makedirs(os.path.dirname(local_path), exist_ok=True)
        try:
            with open(local_path, 'wb') as f:
                self.ftp.retrbinary("RETR {0}".format(ftp_name), f.write)
            # print("downloaded: {0}".format(local_path))
            return True
        except FileNotFoundError:
            print("FAILED FNF: {0}".format(local_path))
        except EOFError:
            print("FAILED EOF: {0}".format(local_path))
            time.sleep(1)
            self.reconnect()
            print("retrying...")
            r = self.download_ftp_file(ftp_name, local_path)
            print("downloaded on retry: {0}".format(local_path))
            return r
        except BrokenPipeError:
            print("FAILED BP: {0}".format(local_path))
        except IsADirectoryError:
            ## same as 550 but pre-exists
            ## os.remove(local_path)
            return 550
        except ftplib.error_perm:
            ## ftplib.error_perm: 550 genomes/TOOLS/EGAP/ortholog_references/9606/current: Not a regular file
            ## its a symlink to a dir.
            os.remove(local_path)
            return 550
        return False

    # item: ('Eublepharis_macularius', {'modify': '20240125225739', 'perm': 'fle', 'size': '4096', 'type': 'dir', 'unique': '6CU599079F', 'unix.group': '562', 'unix.mode': '0444', 'unix.owner': '14'}
    def should_download_file(self, ftp_item, local_name):
        metadata = ftp_item[1]
        ftp_modify = datetime.datetime.strptime(metadata['modify'], '%Y%m%d%H%M%S')
        ftp_size = int(metadata['size'])
        ftp_type = metadata['type']

        local_stat = []
        try:
            local_stat = os.stat(local_name)
        except FileNotFoundError:
            return True
        except NotADirectoryError:
            return True
        # local_is_dir = stat.S_ISDIR(local_stat.st_mode)
        #print(f"item: {ftp_item}")
        #print(f"stat: {local_stat}") 
 
        local_stat_dt = datetime.datetime.fromtimestamp(local_stat.st_mtime)

        #print(f"should_dl: {ftp_size != local_stat.st_size}  {ftp_modify > local_stat_dt}  ")

        if (ftp_type == 'file' and ftp_size != local_stat.st_size) or (ftp_type=='OS.unix=symlink' and ftp_size >= local_stat.st_size):
            return True

        if ftp_modify > local_stat_dt:
            return True

        return False


    ## types
    ## cdir and pdir: . and ..: do nothing
    ## file: a file
    ## dir: a dir
    ## OS.unix=symlink: symlink to a file, treat like a file
    def download_ftp_dir(self, ftp_path, local_path):
        """ replicates a directory on an ftp server recursively """
        for ftp_item in self.ftp.mlsd(ftp_path):
            # print(f"ftp_item: {ftp_item}")
            ##print(f"  name: {ftp_item[0]}")
            ##print(f"  type: {ftp_item[1]['type']}")
            name = ftp_item[0]
            next_ftp_name="/".join([ftp_path,name])
            next_local_name=os.sep.join([local_path,name])
            ftp_type = ftp_item[1]['type']
            ##ftp_types.add(ftp_type)
            if ftp_type=='dir':
                self.download_ftp_dir(next_ftp_name, next_local_name)
            elif ftp_type=='file' or ftp_type=='OS.unix=symlink':
                if self.should_download_file(ftp_item, next_local_name):
                    r = self.download_ftp_file(next_ftp_name, next_local_name)
                    if r == 550:
                        self.download_ftp_dir(next_ftp_name, next_local_name)
                    ##time.sleep(0.1)
            ## else: . or .. do nothing

    def list_ftp_dir(self, ftp_path):
        fl = []
        try: 
            fl = self.ftp.mlsd(ftp_path)
        except:
            return list()
        rl = list()
        for ftp_item in fl:
            ftp_type = ftp_item[1]['type']
            if ftp_type=='dir' or ftp_type=='cdir' or ftp_type=='pdir':
                continue
            name = ftp_item[0]
            rl.append(name)
        return rl


def download_egapx_ftp_data(cache_dir):
    manifest_url = f"{FTP_EGAP_ROOT}/{DATA_VERSION}.mft"
    manifest = urlopen(manifest_url)
    manifest_path = os.path.join(cache_dir, DATA_VERSION + ".mft")
    manifest_list = []
    for line in manifest:
        line = line.decode("utf-8").strip()
        if not line or line[0] == '#':
            continue
        manifest_list.append(line)
        print(f"Downloading {line}")
        ftpd = FtpDownloader()
        ftpd.connect(FTP_EGAP_SERVER)
        ftpd.download_ftp_dir(f"{FTP_EGAP_ROOT_PATH}/{line}", os.path.join(cache_dir, line))
    if cache_dir:
        with open(manifest_path, 'wt') as f:
            for line in manifest_list:
                f.write(f"{line}\n")
    return 0


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


def convert_value(value, key, strict):
    "Convert paths to absolute paths when possible, look into objects and lists"
    if isinstance(value, dict):
        return {k: convert_value(v, key, strict) for k, v in value.items()}
    elif isinstance(value, list):
        return [convert_value(v, key, strict) for v in value]
    else:
        if not isinstance(value, str) or value == '' or re.match(r'[a-z0-9]{2,5}://', value):
            # do not convert - this is a URL or empty string or non-file
            return value
        if not os.path.exists(value):
            if strict:
                raise OSError(ENOENT, f"File for parameter '{key}' doesn't exist", value)
            else:
                return value
        # convert to absolute path
        return os.path.abspath(value)


path_inputs = {'genome', 'hmm', 'softmask', 'reads_metadata', 'short_reads_metadata', 'long_reads_metadata', 'organelles',
               'proteins', 'proteins_trusted', 'reads', 'short_reads', 'long_reads',
               'rnaseq_alignments', 'protein_alignments', 'ortho', 'reference_sets', 'prot_denylist',
               'name_cleanup_rules_file', 'gnomon_filtering_scores_file', 'busco_lineage_download'}
def convert_paths(run_inputs):
    "Convert paths to absolute paths where appropriate"
    input_root = run_inputs['input']
    for key in input_root:
        if key in path_inputs:
            strict = not key.endswith('reads')
            input_root[key] = convert_value(input_root[key], key, strict)
    if 'output' in run_inputs:
        run_inputs['output'] = convert_value(run_inputs['output'], 'output', False)


##def flatten_list(nested_list):
##    return [item for sublist in nested_list if sublist is not None for item in sublist  if item is not None ]
def flatten_list(nested_list):
    flattened_data = []
    for element in nested_list:
        if isinstance(element, list):
            flattened_data.extend(flatten_list(element))
        else:
            flattened_data.append(element)
    return flattened_data

def prepare_reads(run_inputs, force=False, fake_metadata=True):
    """Reformat reads input to be in 'fromFilePairs' format expected by egapx, i.e. [sample_id, [read1, read2]]
    Generate reads metadata file with minimal information - paired/unpaired and valid for existing libraries"""
    res, msg, sras, nosra = True, '', list(), list()
    res1, msg1, sra1, nosra1 = prepare_reads_by_type(run_inputs, 'reads', 'short_reads', fake_metadata)
    res, msg, sras, nosra = res & res1, msg + msg1, sras + sra1, nosra + nosra1 
    res1, msg1, sra1, nosra1 = prepare_reads_by_type(run_inputs, 'short_reads', None, fake_metadata)
    res, msg, sras, nosra = res & res1, msg + msg1, sras + sra1, nosra + nosra1
    res1, msg1, sra1, nosra1  = prepare_reads_by_type(run_inputs, 'long_reads', None, fake_metadata)
    res, msg, sras, nosra = res & res1, msg + msg1, sras + sra1, nosra + nosra1 

    ##print(f'sras 367:   {sras}')
    ##print(f'nosra 368:  {nosra}')
    flt = flatten_list(sras)
    sras = list(set(flatten_list(flt)))  ##sras)))
    ##print(f'sras 371:   {sras}') 

    if 'output' in run_inputs and run_inputs['output'] and os.path.exists(run_inputs['output']):    
        outfile_path = Path(run_inputs['output'], 'sra_metadata.dat')
        save_sra_metadata_file(outfile_path, sras, nosra)
        run_inputs['input']['short_reads_metadata'] = str(outfile_path)
        run_inputs['input']['long_reads_metadata'] = str(outfile_path)

    if force:    
        return res, msg, sras, nosra
    ###for key in ['short_reads_query', 'long_reads_query']:
    # Generate report on the reads, check that the number of SRA runs generated by queries is not too large
    for key in ['short_reads', 'long_reads']:
        if key in run_inputs['input']:
            reads = run_inputs['input'][key]
            paired = 0
            unpaired = 0
            for run in reads:
                _, run_files = run
                if len(run_files) == 1:
                    unpaired += 1
                else:
                    paired += 1
            msg += f"Detected {unpaired} single-end runs, {paired} paired-end runs for query '{key}'\n"
        else:
            query = ''
            # Is it actual Entrez query or synthesized from SRA accession list
            real_query = False
            if key + '_ids' in run_inputs['input']:
                reads = run_inputs['input'][key + '_ids']
                query = "[Accession] OR ".join(reads) + "[Accession]"
                del run_inputs['input'][key + '_ids']
                run_inputs['input'][key + '_query'] = query
                # Should we limit the number of runs similarly to query result?
            elif key + '_query' in run_inputs['input']:
                query = run_inputs['input'][key + '_query']
                real_query = True
            if query:
                sra_uids = sra_uids_query(query)
                sra_runs = sra_metadata_query(sra_uids)
                paired = 0
                unpaired = 0
                #for k, v in sra_runs.items():
                for v in sra_runs:
                    if v[1] == 'paired':
                        paired += 1
                    else:
                        unpaired += 1
                msg += f"Detected {unpaired} single-end runs, {paired} paired-end runs for query '{key}'\n"
                # For actual query check that the returned result is in in limits
                if real_query and len(sra_runs) > RNASEQ_QUERY_LIMIT:
                    if force:
                        msg += 'WARNING: '
                    else:
                        res = False
                        msg += 'ERROR: '
                    msg += f"Number of SRA runs generated by query '{key}' is too large - {len(sra_runs)} with limit {RNASEQ_QUERY_LIMIT}.\n"
                    if not force:
                        msg += "  To ignore this warning and run anyway, use --force option\n"
                
                #sra_uids = sra_uids_query(query)
                #sra_runs = sra_metadata_query(sra_uids_list)
                ##sra_run_ids = sra_query(query)
 
    return res, msg, sras, nosra


def prepare_reads_by_type(run_inputs, reads_type, reads_type_write=None, fake_metadata=True):
    cache_dir = get_cache_dir()
    ##print('cache dir: ', cache_dir)
    res = True
    msg = ''
    sra_to_file_map = {}
    query_to_accessions_map = {}
    sra_metadata_map = {}
    if cache_dir:
        sra_runs_file = os.path.join(cache_dir, SRA_DOWNLOAD_FOLDER, SRA_RUNS_FILE)
        if(os.path.isfile(sra_runs_file)):
            sra_to_file_map = yaml.safe_load(open(sra_runs_file, 'r'))
        sra_queries_file = os.path.join(cache_dir, SRA_DOWNLOAD_FOLDER, SRA_QUERIES_FILE)
        if(os.path.isfile(sra_queries_file)):
            query_to_accessions_map = yaml.safe_load(open(sra_queries_file, 'r'))
        sra_metadata_file = os.path.join(cache_dir, SRA_DOWNLOAD_FOLDER, SRA_METADATA_FILE)
        if(os.path.isfile(sra_metadata_file)):
            sra_metadata_map = yaml.safe_load(open(sra_metadata_file, 'r'))
    #print('sfm: ', sra_to_file_map)
    #print('qam:', query_to_accessions_map)
    #print('sam:', sra_metadata_map)
    
    if reads_type not in run_inputs['input']:
        # print('ln ', 464, reads_type, (reads_type not in run_inputs['input']), ('output' not in run_inputs)  )
        return res, msg, list(), list()
    if reads_type_write is None:
        reads_type_write = reads_type

    prefixes = defaultdict(list)
    sra_to_return = list()
    has_files = False
    reads = run_inputs['input'][reads_type]
    del run_inputs['input'][reads_type]
    ##print('reads_type: ', reads_type)
    ##print("reads: ", reads, " ", str(type(reads)))
    if type(reads) == str:
        # 'reads' is a file name or query
        if os.path.exists(reads):
            # 'reads' is a file name of a tab delimited file of sample names and read files
            # e.g. sample_name read1 read2
            # The file names are either absolute or relative to the current directory
            has_files = True
            with open(reads) as f:
                # Parse lines into run name and read file
                for line in f:
                    line = line.strip()
                    if not line or line[0] == '#':
                        continue
                    mo = re.match(r'([^ \t]+)[ \t]+(.*)', line)
                    if mo:
                        sra_to_return.append(mo.group(2))
                        prefixes[mo.group(1)].append(mo.group(2))
        else:
            # 'reads' is a query
            hash = hashlib.sha1(reads.encode()).hexdigest()
            if query_to_accessions_map and hash in query_to_accessions_map:
                accessions = query_to_accessions_map.get(hash, {})
                if accessions:
                    runs = accessions.get('runs', [])
                    for run in runs:
                        fnames = sra_to_file_map.get(run)
                        if fnames:
                            has_files = True
                            prefixes[run] = fnames
            else:
                uids = sra_uids_query(reads)
                metadata = sra_metadata_query(uids)
                if metadata:
                    for m in metadata:
                        sra_to_return.append(m[0])
                run_inputs['input'][f'{reads_type_write}_query'] = reads
                return res, msg, sra_to_return, list()
    else:
        # 'reads' is a list of files or list of lists
        for rf in reads:
            ##print(f'rf 513 : {rf} :  {len(rf)} : {str(type(rf))} ')
            if type(rf) == str:
                # rf is a file or SRA name
                name = Path(rf).parts[-1]
                mo = re.match(r'([^._]+)', name)
                if mo and mo.group(1) != name:
                    # file name has a prefix separated by '.' or '_', this prefix is the sample name
                    has_files = True
                    sra_to_return.append(mo.group(1))
                    prefixes[mo.group(1)].append(rf)
                else:
                    # TODO: add check that len(Path(rf).parts) == 1, otherwise it is unrecognized path
                    # rf is a SRA name
                    sra_to_return.append(rf)
                    fnames = sra_to_file_map.get(rf)
                    if fnames:
                        has_files = True
                        prefixes[rf] = fnames
            elif type(rf) == list:
                has_files = True
                
                ##print(f'rf 534 : {rf} :  {len(rf)} : {str(type(rf))} ')
                ##print(f'rf 534 : {rf[0]} :  {len(rf[0])} : {str(type(rf[0]))} ')
                ##print(f'rf 534 : {rf[1]} :  {len(rf[1])} : {str(type(rf[1]))} ')
                # find if rf is a tuple already in 'fromFilePairs' format, that is [sample_id, [read1, read2]]
                if len(rf) == 2 and type(rf[0]) == str and type(rf[1]) == list:
                    run_name, run_files = rf
                    ##print(f'540 prefixes : {prefixes}')
                    # Check that run name not yet used
                    if run_name in prefixes:
                        res = False
                        msg += f"Run name {run_name} used multiple times\n"
                        return res, msg, sra_to_return, None
                    else:
                        uids = sra_uids_query(run_name)
                        metadata = sra_metadata_query(uids)
                        ##print(f'549 md : {metadata}')
                        if run_name.startswith('SRR'): 
                            sra_to_return.append(run_name)
                        prefixes[run_name] = list()
                    for r in run_files:
                        if type(r) != str:
                            res = False
                            msg += f"Invalid read input {r}\n"
                            return res, msg, sra_to_return, list()
                        prefixes[run_name].append(r)
                else:
                    # rf is a list of files, we follow star_wnode algorithm to find 'sample_id'
                    names = list(map(lambda x: re.match(r'([^.]+)', Path(x).parts[-1]).group(1), rf))
                    # Find common prefix
                    names.sort()
                    if len(names) == 1:
                        sample_id = names[0]
                    else:
                        for i in range(min(len(names[0]), len(names[-1]))):
                            if names[0][i] != names[-1][i]:
                                break
                        sample_id = names[0][0:i]
                    # Don't end prefix with . or _
                    i = len(sample_id) - 1
                    while i > 0 and sample_id[-1] in "._":
                        sample_id = sample_id[:-1]
                        i -= 1
                    prefixes[sample_id] = rf
    
    ###return res, msg, sra_to_return
    # Create fake metadata file for reads containing only one meaningful information piece - pairedness
    # Have an entry in this file only for SRA runs - files with prefixes SRR, ERR, and DRR and digits
    # Non-SRA reads are handled by rnaseq_collapse internally
    not_sra = list()
    if has_files:
        # run_inputs['input'][f'{reads_type_write}_metadata'] = f.name
        # Always create metadata file even if it's empty
        if fake_metadata: 
            #with tempfile.NamedTemporaryFile(mode='w', delete=False, dir=output, prefix='egapx_reads_metadata_', suffix='.tsv') as f:
            for k, v in prefixes.items():
                # Check file names in v for matching SRR pattern
                mos = [ re.search(r'([SDE]RR[0-9]+)', x) for x in v ]
                # Distinct SRRs
                srrs = list(set([ mo.group(1) for mo in mos if mo ]))
                if len(srrs) > 1:
                    print(f"Error: Run set {k} has multiple run ids in read files: {srrs}")
                    sys.exit(1)
                if srrs:
                    k = srrs[0]
                if re.fullmatch(r'([DES]RR[0-9]+)', k):
                    paired = 'paired' if len(v) == 2 else 'unpaired'
                    # SRR9005248	NA	paired	2	2	NA	NA	NA	NA	NA	NA	NA	0
                    if k in sra_metadata_map:
                        out_rec = list(sra_metadata_map[k]) 
                    else:   
                        out_rec = [k, 'NA', paired, '2', '2', 'NA', 'NA', 'SAMN_'+k, 'NA', 'NA', 'NA', 'NA', '0']
                    #out_rec = [k, 'NA', paired, '2', '2', 'NA', 'NA', , 'NA', 'NA', 'NA', 'NA', '0']
                    not_sra.append(out_rec)
                    # rec = "\t".join(out_rec)
                    ##f.write(rec + '\n')
            ##    f.flush()
            ##    run_inputs['input'][f'{reads_type_write}_metadata'] = f.name
            for k, v in prefixes.items():
                if len(v) > 2:
                    res = False
                    msg += f"More than 2 files for run {k}\n"
                    return res, msg, sra_to_return, not_sra
            run_inputs['input'][reads_type_write] = [ [k, v] for k, v in prefixes.items() ]

    elif reads:
        run_inputs['input'][f'{reads_type_write}_query'] = "[Accession] OR ".join(reads) + "[Accession]"
    
    return res, msg, sra_to_return, not_sra

def get_symbol_format_class_for_taxid(taxid):
    #print('e:get_symbol_format_class_for_taxid')
    #print(f'taxid: {taxid}')
    lineage = get_lineage(taxid)
    #print(f'lineage: {lineage}')

    format_class = ''

    #test_taxids = [7898, 9989, 7742, 6656, 33090]
    #for tt in test_taxids:
    #    print(f'  ct:  {tt} {lineage.count(tt)}')
   
    is_fish = [i for i in lineage if i in FISH]
    is_invert = (ANIMALS in lineage) and not (VERTEBRATE in lineage)
    
    #print(f'is_fish: {is_fish}')
    #print(f'is_invert: {is_invert}')
    #print(f'ANIM: {(ANIMALS in lineage)}')
    #print(f'VERT: {(VERTEBRATE in lineage)}')
    #print(f'INSE: {(INSECTA in lineage)}')
    #print(f'ARTH: {(ARTHROPOD in lineage)}')
    #print(f'LEPI: {(LEPIDOSAURS in lineage)}')
    #print(f'AMPH: {(AMPHIBIANS in lineage)}')

    format_class = 'allupper'
    if is_invert and (INSECTA in lineage or ARTHROPOD in lineage):
        format_class = 'NULL'
    elif (RODENT in lineage) or is_invert:
        format_class = 'uplow'
    elif (LEPIDOSAURS in lineage) or (AMPHIBIANS in lineage) or is_fish:
        format_class = 'alllow'

    #print('x:get_symbol_format_class_for_taxid')
    return format_class


def expand_and_validate_params(run_inputs):
    """ Expand implicit parameters and validate inputs
    Args:
        run_inputs (dict): Input parameters
    Returns:
        bool: True if parameters are valid, False otherwise
    """
    inputs = run_inputs['input']

    if 'taxid' not in inputs:
        print("ERROR: Missing parameter: 'taxid'")
        return False

    taxid = inputs['taxid']
    
    supported, msg = check_supported_taxid(taxid)
    if not supported: 
        print(msg)
        return False
    
    if 'symbol_format_class' not in inputs:
        symbol_format_class = get_symbol_format_class_for_taxid(taxid)
        inputs['symbol_format_class'] = symbol_format_class

    if 'name_cleanup_rules_file' not in inputs:
        name_cleanup_rules_file = get_file_path('misc', 'name_cleanup_rules.txt') 
        inputs['name_cleanup_rules_file'] = name_cleanup_rules_file
    
    if 'gnomon_filtering_scores_file' not in inputs:
        scores_files = list()
        scores_files.append( get_file_path('misc', 'added_scores_general.txt') )
        scores_files.append( get_file_path('misc', 'added_scores_RNAseq.txt') )
        inputs['gnomon_filtering_scores_file'] = scores_files

    if 'genome' not in inputs:
        print("ERROR: Missing parameter: 'genome'")
        return False
    
      
    # Check for proteins input and if empty or no input at all, add closest protein bag
    if 'proteins' not in inputs:
        proteins,trusted = get_closest_protein_bag(taxid)
        if not proteins:
            # Proteins are not specified explicitly and not found by taxid
            print(f"ERROR: Proteins are not found for tax id {inputs['taxid']}")
            print("  This organism is not supported by current NCBI protein database")
            print("  You can specify proteins manually using 'proteins' parameter")
            return False
        inputs['proteins'] = proteins
        if trusted:
            inputs['proteins_trusted'] = trusted

    short_read_keys = {'short_reads', 'short_reads_ids', 'short_reads_query',
                       'reads', 'reads_ids', 'reads_query'}
    short_reads_defs = inputs.keys() & short_read_keys
    if not short_read_keys:
        if inputs['proteins']:
            print("WARNING: It is strongly advised to specify RNA-seq reads using 'short_reads' parameter\n")
        else:
            print("ERROR: Either proteins or RNA-seq reads should be provided for annotation")
            return False
    elif len(short_reads_defs) > 1:
        print("ERROR: Multiple RNA-seq reads definitions are not supported")
        print(f"Remove one of inputs: {', '.join(short_reads_defs)}")
        return False

    train_hmm = False
    if 'hmm' not in inputs:
        best_hmm, good_match = get_closest_hmm(taxid)
        inputs['hmm'] = best_hmm
        train_hmm = not good_match

    if 'train_hmm' not in inputs:
        inputs['train_hmm'] = train_hmm

    if 'max_intron' not in inputs:
        max_intron, genome_size_threshold = get_max_intron(taxid)
        inputs['max_intron'] = max_intron
        inputs['genome_size_threshold'] = genome_size_threshold
    else:
        # Given max_intron is a hard limit, no further calculation is necessary
        inputs['genome_size_threshold'] = 0

    if 'ortho' not in inputs or inputs['ortho'] is None or len(inputs['ortho']) < 4:
        ortho_files = dict()
        if 'ortho' in inputs and isinstance(inputs['ortho'], dict):
            ortho_files = inputs['ortho']

        chosen_taxid=0
        if 'taxid' in ortho_files:
            chosen_taxid = ortho_files['taxid']
        if chosen_taxid == 0: 
            chosen_taxid = get_closest_ortho_ref_taxid(taxid)
        ortho_files['taxid'] = chosen_taxid

        file_id = ['genomic.fna', 'genomic.gff', 'protein.faa']
        for fi in file_id:
            ortho_files[fi] = get_file_path('ortholog_references', f'{chosen_taxid}/current/{fi}.gz')

        ortho_files['name_from.rpt'] = get_file_path('ortholog_references',f'{chosen_taxid}/name_from_ortholog.rpt')
        inputs['ortho'] = ortho_files

    if 'reference_sets' not in inputs or inputs['reference_sets'] is None:
        inputs['reference_sets'] = get_file_path('reference_sets', 'swissprot.asnb.gz')
        inputs['prot_denylist'] = get_file_path('reference_sets', 'swissprot_organelle_bacteria.gi')

    if 'cmsearch' not in inputs or inputs['cmsearch'] is None:
        inputs['cmsearch'] = {'files': [get_file_path('cmsearch', f) for f in "Rfam.seed rfam1410.cm rfam1410_amendments.xml".split()] }
    assert len(inputs['cmsearch']['files']) == 3

    return True


def get_workdir(args):
    workdir = ''
    workdir_file = f"work_dir_{args.executor}.last"
    if args.workdir:
        workdir = args.workdir
        with open(workdir_file, 'w') as f:
            f.write(args.workdir)
    else:
        if os.path.exists(workdir_file):
            with open(workdir_file) as f:
                workdir = f.read().strip()
        else:
            # For cloud-based execution the workdir should be set
            if args.executor in {'aws', 'gcp', 'azure'}:
                print("Work directory not set, use -w at least once")
                return ''
            else:
                workdir = 'work'
    return workdir


def get_cache_dir():
    global user_cache_dir
    if user_cache_dir and os.path.exists(user_cache_dir):
        return user_cache_dir
    return ""


data_version_cache = {}
def get_versioned_path(subsystem, filename):
    """
    Retrieve the versioned path for a given subsystem and filename.

    Checks the cache for the subsystem version, and if not available,
    downloads the manifest file to populate the cache. Returns the path
    including the version if found; otherwise, returns the default path.
    """
    global data_version_cache
    cache_dir = get_cache_dir()
    if not data_version_cache:
        manifest_path = f"{cache_dir}/{DATA_VERSION}.mft"
        if cache_dir and os.path.exists(manifest_path):
            with open(manifest_path, 'rt') as f:
                for line in f:
                    line = line.strip()
                    if not line or line[0] == '#':
                        continue
                    parts = line.split('/')
                    if len(parts) == 2:
                        data_version_cache[parts[0]] = parts[1]
        else:
            manifest_url = f"{FTP_EGAP_ROOT}/{DATA_VERSION}.mft"
            manifest = urlopen(manifest_url)
            manifest_list = []
            for line in manifest:
                line = line.decode("utf-8").strip()
                if not line or line[0] == '#':
                    continue
                parts = line.split('/')
                if len(parts) == 2:
                    data_version_cache[parts[0]] = parts[1]
                    manifest_list.append(line)
            if cache_dir:
                with open(manifest_path, 'wt') as f:
                    for line in manifest_list:
                        f.write(f"{line}\n")

    if subsystem not in data_version_cache:
        return os.path.join(subsystem, filename)
    version = data_version_cache[subsystem]
    return os.path.join(subsystem, version, filename)


# Get file path for subsystem, either from cache or from remote server
def get_file_path(subsystem, filename):
    cache_dir = get_cache_dir()
    vfn = get_versioned_path(subsystem, filename)
    file_path = os.path.join(cache_dir, vfn)
    file_url = f"{FTP_EGAP_ROOT}/{vfn}"
    if os.path.exists(file_path):
        return file_path
    return file_url


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


lineage_cache = {}
taxname_cache = {}
def get_lineage_with_ranks(taxid):
    global lineage_cache
    global taxname_cache
    ranks = {}
    if not taxid:
        return [], {}
    if taxid in lineage_cache:
        return lineage_cache[taxid]
    # Try cached taxonomy database
    taxonomy_db_name = os.path.join(get_cache_dir(), get_versioned_path("taxonomy", "taxonomy.sqlite3"))
    if os.path.exists(taxonomy_db_name):
        conn = sqlite3.connect(taxonomy_db_name)
        if conn:
            c = conn.cursor()
            lineage = [taxid]
            cur_taxid = taxid
            while cur_taxid != 1:
                c.execute("SELECT parent, rank FROM TaxidInfo WHERE taxid = ?", (cur_taxid,))
                parent, rank = c.fetchone()
                if rank:
                    ranks[cur_taxid] = rank
                lineage.append(parent)
                cur_taxid = parent
            lineage.reverse()
            for t in lineage:
                if t not in taxname_cache:
                    c.execute("SELECT name FROM TaxNames WHERE taxid = ?", (t,))
                    name = c.fetchone()[0]
                    taxname_cache[t] = name
            c.close()
            conn.close()
            lineage_cache[taxid] = lineage, ranks
            return lineage, ranks
    
    # Fallback to API
    taxon_json_file = urlopen(dataset_taxonomy_url+str(taxid))
    taxon = json.load(taxon_json_file)["taxonomy_nodes"][0]
    # print(f"taxid: {taxid}, taxon: {taxon}")
    lineage = taxon["taxonomy"].get("lineage", [])
    lineage.append(taxon["taxonomy"]["tax_id"])

    ranks_file = urlopen(dataset_taxonomy_url+",".join(map(str, lineage)))
    nodes = json.load(ranks_file)["taxonomy_nodes"]
    for node in nodes:
        if 'rank' in node["taxonomy"]:
            ranks[node["taxonomy"]["tax_id"]] = node["taxonomy"]["rank"]
        if 'organism_name' in node["taxonomy"]:
            taxname_cache[node["taxonomy"]["tax_id"]] = node["taxonomy"]["organism_name"]

    lineage_cache[taxid] = lineage, ranks
    return lineage, ranks


def get_lineage(taxid):
    lineage, _ = get_lineage_with_ranks(taxid)
    return lineage


def get_tax_name(taxid):
    """
    Returns the taxonomic name for a given taxid.
    NB: the name cache is valid only AFTER getting the lineage for the taxid or one of its descendants
    """
    global taxname_cache
    return taxname_cache.get(taxid)


def get_tax_file(subsystem, tax_path):
    vfn = get_versioned_path(subsystem, tax_path)
    taxids_path = os.path.join(get_cache_dir(), vfn)
    taxids_url = f"{FTP_EGAP_ROOT}/{vfn}"
    taxids_file = []
    if os.path.exists(taxids_path):
        with open(taxids_path, "rb") as r:
            taxids_file = r.readlines()
    else:
        taxids_file = urlopen(taxids_url)
    return taxids_file


def check_supported_taxid(taxid):
    if not taxid:
        return False, "No taxid provided"

    lineage = get_lineage(taxid)
    reqs = [ARTHROPOD, VERTEBRATE, MAGNOLIOPSIDA, ECHINODERMATA]
    lineage_check = any(l in lineage for l in reqs)

    #print(f'input tax: {taxid} , lineage: {lineage} , contains: { lineage_check } ')
    if lineage_check:
        return True, ''
    return False, "ERROR: Input taxid not within supported range. Must be under Arthropoda(6656), Vertebrata(7742), Magnoliopsida(3398), or Echinodermata(7586) according to NCBI Taxonomy (https://www.ncbi.nlm.nih.gov/taxonomy)"


def get_closest_protein_bag(taxid):
    if not taxid:
        return '',''

    taxids_file = get_tax_file("target_proteins", "taxid.list")
    taxids_list = []
    for line in taxids_file:
        line = line.decode("utf-8").strip()
        if len(line) == 0 or line[0] == '#':
            continue
        parts = line.split('\t')
        if len(parts) > 0:
            t = parts[0]
            taxids_list.append(int(t))

    lineage = get_lineage(taxid)

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
        return '',''
    # print(best_taxid, best_score)
    protein_file = get_file_path("target_proteins", f"{best_taxid}.faa.gz")
    trusted_file = get_file_path("target_proteins", f"{best_taxid}.trusted_proteins.gi")
    return protein_file,trusted_file


def get_closest_hmm(taxid):
    """
    Given a taxid, this function returns the closest HMM parameters file and a boolean indicating whether the match is good or not.
    
    Parameters:
        taxid (int): The taxid for which to find the closest HMM parameters file.
    
    Returns:
        tuple: A tuple containing two elements:
            - str: The path to the closest HMM parameters file.
            - bool: A boolean indicating whether the match is good or not.
    
    Raises:
        None
    
    Examples:
        >>> get_closest_hmm(123456)
        ('/path/to/hmm_parameters/123456.params', True)
    """
    if not taxid:
        return "", False

    taxids_file = get_tax_file("gnomon", "hmm_parameters/taxid.list")

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

    lineage, ranks = get_lineage_with_ranks(taxid)

    is_mammal = MAMMAL in lineage
    best_lineage = None
    best_taxid = None
    best_score = 0
    best_good_match = False
    for (t, l) in lineages:
        pos1 = 0
        last_match = 0
        last_good_match = False
        for pos in range(len(lineage)):
            tax_id = lineage[pos]
            rank = ranks.get(tax_id, '')
            while tax_id != l[pos1]:
                if pos1 + 1 < len(l):
                    pos1 += 1
                else:
                    break
            if tax_id == l[pos1]:
                last_match = pos1
                last_good_match = last_good_match or rank == 'GENUS' or (rank == 'FAMILY' and is_mammal)
            else:
                break
        if last_match > best_score:
            best_score = last_match
            best_taxid = t
            best_lineage = l
            best_good_match = last_good_match

    if best_score == 0:
        return "", False

    # Either perfect match or match of the genus or family for mammals
    good_match = best_good_match or best_score == len(best_lineage) - 1
    # print(best_lineage)
    # print(best_taxid, best_score)
    return get_file_path("gnomon", f"hmm_parameters/{best_taxid}.params"), good_match


def get_closest_ortho_ref_taxid(taxid):
    if not taxid:
        return 0
 
    taxids_file = get_tax_file("ortholog_references", "taxid.list")

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

    lineage, ranks = get_lineage_with_ranks(taxid)
    
    is_mammal = MAMMAL in lineage
    best_lineage = None
    best_taxid = None
    best_score = 0
    best_good_match = False
    for (t, l) in lineages:
        pos1 = 0
        last_match = 0
        last_good_match = False
        for pos in range(len(lineage)):
            tax_id = lineage[pos]
            rank = ranks.get(tax_id, '')
            while tax_id != l[pos1]:
                if pos1 + 1 < len(l):
                    pos1 += 1
                else:
                    break
            if tax_id == l[pos1]:
                last_match = pos1
                last_good_match = last_good_match or rank == 'GENUS' or (rank == 'FAMILY' and is_mammal)
            else:
                break
        if last_match > best_score:
            best_score = last_match
            best_taxid = t
            best_lineage = l
            best_good_match = last_good_match
            

    if best_score == 0:
        return 0
    
    # Either perfect match or match of the genus or family for mammals
    # good_match = best_good_match or best_score == len(best_lineage) - 1
    ##print(best_lineage)
    ##print(best_taxid, best_score)
    return best_taxid


def download_busco_lineage_set(cache_dir):
    busco_downloads_dir = os.path.join(cache_dir, "busco_downloads")
    busco_list_file = os.path.join(busco_downloads_dir, "file_versions.tsv")
    busco_list_url = f"{BUSCO_DATA_URL}/file_versions.tsv"
    if not os.path.exists(busco_list_file):
        os.makedirs(busco_downloads_dir, exist_ok=True)
        urllib.request.urlretrieve(busco_list_url, busco_list_file)


def download_busco_lineage(lineage, cache_dir):
    busco_downloads_dir = os.path.join(cache_dir, "busco_downloads")
    lineage_path = os.path.join(busco_downloads_dir, "lineages", lineage)
    if os.path.isdir(lineage_path):
        return
    lineage_set = get_busco_lineage_set()
    lineage_date = lineage_set.get(lineage)
    if not lineage_date:
        raise Exception(f"Lineage {lineage} not found")
    # Download the lineage tar.gz and unpack it
    busco_lineage_url = f"{BUSCO_DATA_URL}/lineages/{lineage}.{lineage_date}.tar.gz"
    busco_lineage_file = os.path.join(busco_downloads_dir, f"{lineage}.tar.gz")
    os.makedirs(busco_downloads_dir, exist_ok=True)
    urllib.request.urlretrieve(busco_lineage_url, busco_lineage_file)
    with tarfile.open(busco_lineage_file) as tar:
        tar.extractall(path=os.path.join(busco_downloads_dir, "lineages"))
    os.remove(busco_lineage_file)


def get_busco_lineage_set():
    """
    Returns a set of all BUSCO lineages
    """
    cache_dir = get_cache_dir()
    busco_list_file = os.path.join(cache_dir, "busco_downloads", "file_versions.tsv")
    busco_list_url = f"{BUSCO_DATA_URL}/file_versions.tsv"
    if os.path.exists(busco_list_file):
        return { x[0]: x[1] for x in map(lambda x: x.split('\t'), open(busco_list_file, 'r').readlines())}
    else:
        return { x[0]: x[1] for x in map(lambda x: x.decode("utf-8").split('\t'), urlopen(busco_list_url).readlines())}


busco_name_fixes = {
    'artiodactyla': 'cetartiodactyla',
    'eudicotyledons': 'eudicots'
}
def get_busco_lineage(taxid):
    lineage = get_lineage(taxid)
    busco_lineage_set = get_busco_lineage_set()
    for t in lineage.__reversed__():
        name = get_tax_name(t).lower()
        name = busco_name_fixes.get(name, name) + '_odb10'
        if name in busco_lineage_set:
            return name
    return None


def get_max_intron(taxid):
    if not taxid:
        return 0, 0
    lineage = get_lineage(taxid)
    if VIRIDIPLANTAE in lineage:
        return 300000, 3000000000
    elif VERTEBRATE in lineage:
        return 1200000, 2000000000
    else:
        return 600000, 500000000


def to_dict(x: List[str]):
    d = {}
    s = len(x)
    i = 0
    while i < s:
        el = x[i]
        i += 1
        if el and el[0] == '-':
            if i < s:
                v = x[i]
                if v and (v[0] != '-'  or (v[0] == '-' and ' ' in v)):
                    d[el] = v
                    i += 1
                else:
                    d[el] = ""
            else:
                d[el] = ""
        else:
            d[el] = ""
    return d

def merge_params(task_params, run_inputs):
    # Walk over the keys in run_inputs and merge them into task_params
    for k in run_inputs.keys():
        if k in task_params:
            if isinstance(task_params[k], dict) and isinstance(run_inputs[k], dict):
                task_params[k] = merge_params(task_params[k], run_inputs[k])
            else:
                task_params_dict = to_dict(shlex.split(task_params[k]))
                run_inputs_dict = to_dict(shlex.split(run_inputs[k]))
                task_params_dict.update(run_inputs_dict)
                task_params_list = []
                for k1, v in task_params_dict.items():
                    task_params_list.append(k1)
                    if v:
                        task_params_list.append(v)
                task_params[k] = shlex.join(task_params_list)
        else:
            task_params[k] = run_inputs[k]
    return task_params


filter_out = {
"genes (other)",
"genes (transcribed pseudo)",
"genes (minor correction)",
"mRNAs (has gaps)",
"mRNAs (model with correction)",
"non-coding RNAs (partial)",
"non-coding RNAs (correction)",
"non-coding RNAs (ab initio > 5%)",
"non-coding RNAs (has gaps)",
"pseudo transcripts (partial)",
"pseudo transcripts (correction)",
"pseudo transcripts (has gaps)",
"CDSs (exon <= 3nt)",
"CDSs (minor correction)",
}
def print_statistics(output):
    feature_counts_txt = Path(output) / 'stats' / 'feature_counts.txt'
    if not feature_counts_txt.exists():
        return 1
    with open(feature_counts_txt, 'rt') as f:
        for line in f:
            line = line.rstrip()
            parts = line.split(': ')
            if len(parts) == 2:
                if parts[0].strip() not in filter_out:
                    print(f"{parts[0]}: {parts[1]}")
            else:
                print(line)
    return 0


def print_statistics_old(output):
    # Old way
    accept_gff = Path(output) / 'complete.genomic.gff'
    print(f"Statistics for {accept_gff}")
    counter = defaultdict(int)
    if accept_gff.exists():
        with open(accept_gff, 'rt') as f:
            for line in f:
                line = line.strip()
                if not line or line[0] == '#':
                    continue
                parts = line.split()
                if len(parts) < 3:
                    continue
                counter[parts[2]] += 1
    keys = list(counter.keys())
    keys.sort()
    for k in keys:
        print(f"{k:12s} {counter[k]}")
    return 0


def print_log_values(outdir, task_key, report_effective_values):
    task_dst = Path(outdir) / "execution_logs" / task_key
    logs = glob.glob(os.path.join(task_dst, ".command.*"))
    for log in logs:
        os.rename(log, os.path.join(task_dst, os.path.basename(log)[1:]))
    if task_key in report_effective_values:
        parameter_names = report_effective_values[task_key]
        effective_value_log = Path(task_dst) / "command.log"
        effective_values = []
        with open(effective_value_log, 'rt') as f:
            for line in f:
                parts = line.strip().split(None, 1)
                if len(parts) > 1 and  parts[0] in parameter_names:
                    effective_values.append((parts[0], parts[1]))
        for parameter_name, effective_value in effective_values:
            print(f"Effective parameter {parameter_name} is {effective_value}")


# Map of tasks to report effective values, key is the task name, value is a set of parameter names
g_report_effective_values = {
    "setup_genome/get_genome_info" : {"max_intron"},
}
def collect_logs(outdir):
    print("Collecting execution logs")
    run_trace_path = Path(outdir) / "run.trace.txt"
    if not os.path.exists(run_trace_path):
        return 1
    with open(Path(outdir) / "run_work_dir.txt", "rt") as f:
        workdir = f.read()
    with open(run_trace_path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line or line[0] == '#':
                continue
            parts = line.split('\t')
            if parts[0] == "task_id":
                continue
            _, hash, _, name = parts[:4]
            name_parts = name.split(":")
            if name_parts[0] == "egapx":
                name_parts = name_parts[1:]
            # if name ends with (number), remove it and add number as a suffix

            # AA: the tag in parentheses could contain an arbitrary string,
            # not just a a number. e.g. "my_task (foo 42)". In such case 
            # task_dst will not be a valid path, containing parentheses and spaces,
            # which will cause the logic below to crash.
            #
            # Instead, the path needs to be cleaned up (e.g. strip garbarge and 
            # replace non-alphanumeric characters with underscores)
            #
            # For now I'm going to adjust the nextflow processes to conform to
            # the expectations here, because I don't know if any downstream logic
            # relies on that as well.

            mo = re.match(r'(.*)\s+\((\d+)\)$', name_parts[-1])
            if mo:
                name_parts[-1] = f"{mo.group(1)}_{mo.group(2)}"
            task_key = os.path.join(*name_parts)
            task_dst = Path(outdir) / "execution_logs" / task_key
            os.makedirs(task_dst, exist_ok=True)
            # The following code is specific to the cloud used. Two cases, AWS and filesystem are implemented
            # now. Potential future support for other cloud providers is possible.
            # Get full name of task workdir from hash
            hash_parts = hash.split("/")
            if workdir.startswith("s3://"):
                # AWS S3
                task_workdir_stub = os.path.join(workdir, hash)
                cp = subprocess.run(["aws", "s3", "ls", task_workdir_stub], check=True, capture_output=True)
                # parse output in form PRE dirname/
                lines = cp.stdout.decode('utf-8').split('\n')
                # print(lines)
                for line in lines:
                    line = line.strip()
                    if line.startswith("PRE "):
                        task_workdir = os.path.join(workdir, hash_parts[0], line[4:])
                        break
                subprocess.run(["aws", "s3", "cp", task_workdir, task_dst, "--recursive", "--exclude", "*", "--include", ".command.*"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            elif workdir.startswith("gs://"):
                # Google Cloud Storage, not yet implemented
                print("Collecting logs from Google Cloud Storage is not yet implemented")
            elif workdir.startswith("az://"):
                # Azure Blob Storage, not yet implemented
                print("Collecting logs from Azure Blob Storage is not yet implemented")
            else:
                # Local filesystem
                task_workdir_stub = os.path.join(workdir, hash_parts[0])
                cp = subprocess.run(["ls", "-1", task_workdir_stub], check=True, capture_output=True)
                lines = cp.stdout.decode('utf-8').split('\n')
                # print(lines)
                for line in lines:
                    if line.startswith(hash_parts[1]):
                        task_workdir = os.path.join(task_workdir_stub, line)
                        break

                cp = subprocess.run(f"cp {task_workdir}/.command.* {task_dst}", shell=True, check=True, capture_output=True)

            logs = glob.glob(os.path.join(task_dst, ".command.*"))
            # Make logs visible
            for log in logs:
                os.rename(log, os.path.join(task_dst, os.path.basename(log)[1:]))
            print_log_values(outdir, task_key, g_report_effective_values)
            # print(f"collected from {hash}, task {name} {task_workdir} to {task_dst}")
    return 0


def get_software_version():
    global software_version
    if not software_version or software_version[0] == '$':
        process_result = subprocess.run(["git", "describe", "--tags"], stdout=subprocess.PIPE, check=False)
        if process_result.returncode == 0:
            version = process_result.stdout.decode('utf-8').strip()
        else:
            version = "unknown"
    else:
        version = software_version
    return version


def sra_uids_query(sra_potential_query):
    time.sleep(1)
    biomol_str = " AND biomol_transcript\[properties\] "
    
    if 'biomol_transcript' not in sra_potential_query:
        sra_potential_query += biomol_str
   
    esearch = urlopen(esearch_url+quote(sra_potential_query.encode('utf-8') ))
    esearch_json = json.load(esearch)
    uids = esearch_json["esearchresult"]["idlist"]
    return uids
    ##uids = ','.join(esearch_json["esearchresult"]["idlist"])
    ##print(esearch_json)


def sra_metadata_query(sra_uids_list):
    uids = ','.join(sra_uids_list)
    runinfo = urlopen(runinfo_url+uids)
    runinfo_csv = runinfo.read().decode("utf-8")
    ##print(runinfo_csv)

    lines = runinfo_csv.split(NL)
    header = lines[0]
    body = lines[1:]
    keypos = {}
    for i, k in enumerate(header.split(",")):
        keypos[k] = i
        ##print(i, " : ", k)
    ##sra_meta_keys = ["SRA run accession", "SRA sample accession", "Run type", "SRA read count", "SRA base count", "Average insert size", "Insert size stdev",
    ##            "Platform type", "Model", "SRA Experiment accession", "SRA Study accession", "Biosample accession", "Bioproject ID", "Bioproject Accession", "Scientific name", "TaxID", "Release date"]
    
    records = {}
    for parts in csv.reader(body, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True):
        if len(parts) == 0:
            continue
        if "Run" not in keypos:
            continue
        run = parts[keypos["Run"]]
        record = [run]
        paired = False
        for k in ["Sample", "LibraryLayout", "spots", "bases", "NA", "NA", "Platform", "Model", "Experiment", "SRAStudy", "BioSample", "ProjectID", "BioProject", "ScientificName", "TaxID", "ReleaseDate"]:
            if k not in keypos:
                record.append("NA")
                continue
            v = parts[keypos[k]]
            if k == "LibraryLayout":
                if parts[keypos["LibraryLayout"]] == "PAIRED":
                    record.append("paired")
                    paired = True
                else:
                    record.append("unpaired")
            elif k == "spots":
                if not paired:
                    record.append(v)
                else:
                    i = int(v)
                    if 'spots_with_mates' in keypos:
                        i += int(parts[keypos['spots_with_mates']])
                    record.append(str(i))
            else:
                record.append(v)
        records[run] = record
        # print(TAB.join([run] + record))
    return records.values()

def save_sra_metadata_file(dest, sras, non_sras):
    time.sleep(1)
    uids=list()
    metadata=list()
    if(sras):
        try:
            uids = sra_uids_query(" OR ".join(sras))
            metadata = sra_metadata_query(uids)
        except:
            print('eutil connection failure')

    sra_meta_keys = ["SRA run accession", "SRA sample accession", "Run type", "SRA read count", "SRA base count", "Average insert size", "Insert size stdev",
                "Platform type", "Model", "SRA Experiment accession", "SRA Study accession", "Biosample accession", "Bioproject ID", "Bioproject Accession", "Scientific name", "TaxID", "Release date"]
  
    seen_sras = set()
    with open(dest, 'wt') as outf:
        print(f"#{TAB.join(sra_meta_keys)}", file=outf)
        for md in metadata:
            print(TAB.join(md), file=outf)
            seen_sras.add(md[0])
            
        for ns in non_sras:
            if ns is not [] and ns[0] not in seen_sras:
                print(TAB.join(ns), file=outf) 
     

def expand_sra_query(query):
    if not query:
        return [], {}
    # get the list of SRA accessions to be downloaded and metadata if available
    if type(query) == str:
        # either actual query or file name
        if os.path.exists(query):
            return [], {}
        else:
            # Runs an actual Entrez query and returns the list of SRA accessions
            sra_uids = sra_uids_query(query)
            sra_runs = sra_metadata_query(sra_uids)

            metadata = dict()
            for sr in sra_runs:
                run_id = sr[0]
                metadata[run_id] = sr

            return list(metadata.keys()), metadata
            ##return sra_runs, {}        
            ##metadata = sra_query(query)
            ##return list(metadata.keys()), metadata
    elif type(query) == list:
        # list of SRA accessions
        return query, {}
    return query, {}


def download_sra_query(query, sra_dir):
    sra_runs_list, metadata = expand_sra_query(query)
 
    if len(sra_runs_list) == 0:
        return 0
    # check which fasta files are already there and should not be redownloaded
    existing_fasta_files = glob.glob(os.path.join(sra_dir, '*.fasta'))
    # get the list of sra to be downloaded
    find_sra_name = lambda x : os.path.splitext(os.path.basename(x))[0].split('_')[0]
    existing_sras = set(map(find_sra_name, existing_fasta_files))
    sra_list = list(filter(lambda x : x not in existing_sras, sra_runs_list))
    if len(sra_list) == 0:
        return 0
    #write the sra_list to be downloaded in file
    Path(sra_dir).mkdir(parents=True, exist_ok=True)
    for sra in sra_list:
        print(f"Downloading {sra}")
        defline = ">gnl|SRA|$ac.$si.$ri"
        # common arguments for fasterq-dump or fastq-dump
        args = [ sra, '--skip-technical', '--split-files', '-O', sra_dir ]
        fasterq_args = ['fasterq-dump'] + args + ['--fasta', '--seq-defline', defline, '--threads', '6']
        fastq_args = ['fastq-dump'] + args + ['--fasta', '0', '--defline-seq', defline]
        # try:
        # run fasterq-dump or fastq-dump
        if shutil.which('fasterq-dump'):
            try:
                subprocess.run(fasterq_args, check=True)
            except subprocess.CalledProcessError:
                print("Error: fasterq-dump failed")
                try:
                    subprocess.run(fastq_args, check=True)
                except subprocess.CalledProcessError:
                    print("Error: fastq-dump failed")
                    return 1
        elif shutil.which('fastq-dump'):
            try:
                subprocess.run(fastq_args, check=True)
            except subprocess.CalledProcessError:
                print("Error: fastq-dump failed")
                return 1
        else:
            print("Error: neither fasterq-dump nor fastq-dump is installed")
            return 1

    # update sra_reads yaml file with the map from SRA accession
    # to list of fasta files downloaded for this accession
    existing_fasta_files = glob.glob(os.path.join(sra_dir, '*.fasta'))

    sra_runs_fn = os.path.join(sra_dir, SRA_RUNS_FILE)
    if os.path.exists(sra_runs_fn):
        sra_runs_dict = yaml.safe_load(open(sra_runs_fn, 'r'))
    else:
        sra_runs_dict = {}
    for sra in sra_list:
        sra_runs_dict[sra] = [os.path.abspath(sra_file) for sra_file in existing_fasta_files if find_sra_name(sra_file) == sra]
    with open(sra_runs_fn, 'w') as f:
        yaml.dump(sra_runs_dict, f)
        f.flush()
    
    if metadata:
        sra_metadata_fn = os.path.join(sra_dir, SRA_METADATA_FILE)
        if os.path.exists(sra_metadata_fn):
            sra_metadata_dict = yaml.safe_load(open(sra_metadata_fn, 'r'))
        else:
            sra_metadata_dict = {}
        sra_metadata_dict.update(metadata)
        with open(sra_metadata_fn, 'w') as f:
            yaml.dump(sra_metadata_dict, f)
            f.flush()
        
        sra_queries_fn = os.path.join(sra_dir, SRA_QUERIES_FILE)
        if os.path.exists(sra_queries_fn):
            sra_queries_dict = yaml.safe_load(open(sra_queries_fn, 'r'))
        else:
            sra_queries_dict = {}
        hash = hashlib.sha1(query.encode()).hexdigest()
        sra_queries_dict[hash] = { "query": query, "runs": sra_runs_list}
        # sra_queries_dict[query] = sra_reads_list
        with open(sra_queries_fn, 'w') as f:
            yaml.dump(sra_queries_dict, f)
            f.flush()
    return 0


def download_offline_data(args):
    """
    Download offline data necessary for EGAPx operation.

    This function handles the downloading of support data, SRA data, and BUSCO lineage data
    to a specified local cache directory. It respects the dry-run mode and supports command line
    inputs for filenames and output paths. If the input filename is provided, it downloads the
    data for SRA runs and BUSCO lineage data for the specific runs and taxonomy,

    Args:
        args: A namespace object containing command line arguments, expected to have:
            - local_cache: The path to the local cache directory.
            - dry_run: A flag to indicate dry-run mode.
            - filename: The name of the file containing user input.
            - output: The output directory path.
            - force: A flag to force certain operations.

    Returns:
        int: Returns 1 if local cache is not set or an error occurs during download,
             otherwise returns 0 upon successful dry-run.
    """

    # print(f"Download only: {args.download_only}")
    if not args.local_cache:
        print("Local cache not set, please use -lc option")
        return 1
    sra_dir = os.path.abspath(os.path.join(args.local_cache, SRA_DOWNLOAD_FOLDER))
    if args.dry_run:
        print(f"Download support data to {args.local_cache}")
        if args.filename:
            print(f"Download SRA to {sra_dir}")
        return 0
    os.makedirs(args.local_cache, exist_ok=True)
    download_egapx_ftp_data(args.local_cache)
    if args.filename:
        run_inputs = repackage_inputs(yaml.safe_load(open(args.filename, 'r')))
        if args.output:
            run_inputs['output'] = args.output
        print(f"Download SRA to {sra_dir}")
        res, reads_msg, sras, nosras = prepare_reads(run_inputs, args.force, False)
        if not res:
            print('return before downloads')
            print(reads_msg, end='')
            return 1
        else:
            inputs = run_inputs['input']
            short_reads_query = inputs.get('short_reads_query', '')
            long_reads_query  = inputs.get('long_reads_query', '')
            download_sra_query(short_reads_query, sra_dir)
            download_sra_query(long_reads_query, sra_dir)
            # Download BUSCO lineage
            download_busco_lineage_set(args.local_cache)
            if 'busco_lineage' in inputs:
                busco_lineage = inputs['busco_lineage']
            else:
                busco_lineage = get_busco_lineage(inputs['taxid'])
            print(f"Downloading BUSCO lineage {busco_lineage}")
            download_busco_lineage(busco_lineage, args.local_cache)
    return 0


def main(argv):
    "Main script for EGAPx"

    # Parse command line
    args = parse_args(argv)

    if args.report_version:
        print(f"EGAPx {get_software_version()}")
        return 0

    #warn user that this is an alpha release
    print("\n!!WARNING!!\nThis is an alpha release with limited features and organism scope to collect initial feedback on execution. Outputs are not yet complete and not intended for production use.\n")
    
    if args.download_only:
        return download_offline_data(args)

    # Command line overrides manifest input
    run_inputs = dict()
    if args.output:
        run_inputs['output'] = args.output
    
    # Separate modes of execution which depend only on output directory
    if not args.dry_run and (not 'output' in run_inputs or not run_inputs['output']):
        print("Output directory not set")
        return 1

    if args.summary_only:
        if not os.path.exists(run_inputs['output']):
            print(f"Output directory {run_inputs['output']} does not exist")
            return 1
        return print_statistics(run_inputs['output'])
    
    if args.logs_only:
        if not os.path.exists(run_inputs['output']):
            print(f"Output directory {run_inputs['output']} does not exist")
            return 1
        return collect_logs(os.path.join(run_inputs['output'], 'nextflow'))

    # Check that input and output set
    if not args.filename:
        print("Input file must be set")
        return 1

    global user_cache_dir
    if args.local_cache:
        # print(f"Local cache: {args.local_cache}")
        if not os.path.exists(args.local_cache):
            print(f"Local cache directory {args.local_cache} does not exist")
            return 1
        user_cache_dir = args.local_cache

    packaged_distro = bool(getattr(sys, '_MEIPASS', ''))
    script_directory = getattr(sys, '_MEIPASS', os.path.abspath(os.path.dirname(__file__)))
    config_file = get_config(script_directory, args)
    if not config_file:
        return 1

    # Check for workdir, set if not set, and manage last used workdir
    workdir = get_workdir(args)
    if not workdir:
        return 1
   
    # Read default task parameters into a dict
    task_params = yaml.safe_load(open(Path(script_directory) / 'assets' / 'default_task_params.yaml', 'r'))
    run_inputs = repackage_inputs(yaml.safe_load(open(args.filename, 'r')))

    # Command line overrides manifest input
    if args.output:
        run_inputs['output'] = args.output

    # Shortcut for following code    
    inputs = run_inputs['input']

    if int(args.ortho_taxid) > 0:
        inputs['ortho'] = {'taxid': int(args.ortho_taxid)}
    
    if not expand_and_validate_params(run_inputs):
        return 1

    # In GNOMON's chainer module, the default is -minlen 165 and -minscor 25.0,
    # use -minlen 225 and -minscor 40.0 for Magnoliopsida and Vertebrates,
    lineage = get_lineage(inputs['taxid'])
    if MAGNOLIOPSIDA in lineage or VERTEBRATE in lineage:
        minlen = 225
        minscor = 40.0
    else:
        minlen = 165
        minscor = 25.0
    task_params = merge_params(task_params, {'tasks': { 'chainer_wnode': {'chainer_wnode': f"-minlen {minlen} -minscor {minscor}"}}})

    # Add some parameters to specific tasks
    final_asn_params = f"-annot-software-version {get_software_version()}"
    if 'annotation_provider' in inputs:
        final_asn_params += f" -annot-provider \"{inputs['annotation_provider']}\""
    run_date_std = time.strftime("%Y_%m_%d", time.localtime(start_time))
    if 'annotation_name_prefix' in inputs:
        annotation_name_prefix = re.sub(r'\s+', '-', inputs['annotation_name_prefix'])
        final_asn_params += f" -annot-name {annotation_name_prefix}-GB_{run_date_std}"
    else:
        final_asn_params += f" -annot-name GB_{run_date_std}"
    # NB: this date format is required by NCBI
    run_date_us = time.strftime("%m/%d/%Y", time.localtime(start_time))
    final_asn_params += f" -annot-date {run_date_us}"
    if 'locus_tag_prefix' in inputs:
        locus_tag_prefix = re.sub(r'\s+', '-', inputs['locus_tag_prefix'])
        final_asn_params += f" -locus-tag-prefix {locus_tag_prefix}"
    task_params = merge_params(task_params, {'tasks': { 'final_asn_markup': {'final_asn': final_asn_params}}} )

    if 'busco_lineage' in inputs:
        busco_lineage = inputs['busco_lineage']
    else:
        busco_lineage = get_busco_lineage(inputs['taxid'])
        inputs['busco_lineage'] = busco_lineage
    if not busco_lineage:
        # This should never happen, BUSCO lineages cover all taxids
        print(f"BUSCO lineage not found for taxid {inputs['taxid']}")
        return 1
    
    cache_dir = get_cache_dir()
    if cache_dir:
        busco_lineage_file_name = os.path.join(cache_dir, 'busco_downloads', 'lineages', busco_lineage)
        if os.path.exists(busco_lineage_file_name):
            inputs['busco_lineage_download'] = busco_lineage_file_name
        else:
            print(f"BUSCO lineage {busco_lineage} not found in {cache_dir}")
            print(f"Please run egapx.py -lc {cache_dir} -dl {args.filename} to download the lineage")
            return 1

    # Create output directory if needed
    if 'output' in run_inputs and run_inputs['output']: 
        os.makedirs(os.path.join(run_inputs['output'], 'nextflow'), exist_ok=True)

    # Reformat reads into pairs in fromPairs format and add reads_metadata.tsv file
    res, reads_msg, sras, nosras = prepare_reads(run_inputs, args.force)
    if not res:
        print(reads_msg, end='')
        return 1

    # Convert inputs to absolute paths
    try:
        convert_paths(run_inputs)
    except OSError as e:
        print(F"{e.strerror}: {e.filename}")
        return 1

    # Add to default task parameters, if input file has some task parameters they will override the default
    task_params = merge_params(task_params, run_inputs)
    ##exit(1)

    # Move output from YAML file to arguments to have more transparent Nextflow log
    if 'output' in task_params:
        output = task_params['output']
        del task_params['output']
    elif not args.dry_run:
        # Should never happen
        print("Impossible condition: no output directory")
        return 1
    else:
        output = ''
    nextflow_out_dir = os.path.join(output, 'nextflow')

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
    top_level_script = os.environ.get("EGAPX_TOP_LEVEL_SCRIPT", 'ui.nf')
    if packaged_distro:
        main_nf = Path(script_directory) / 'nf' / top_level_script
    else:
        main_nf = Path(script_directory) / '..' / 'nf' / top_level_script
    nf_cmd = ["nextflow", "-C", config_file, "-log", f"{nextflow_out_dir}/nextflow.log", "run", main_nf, "--output", output]
    if args.stub_run:
        nf_cmd += ["-stub-run", "-profile", "stubrun"]
    if args.report:
        nf_cmd += ["-with-report", f"{args.report}.report.html", "-with-timeline", f"{args.report}.timeline.html"]
    else:
        nf_cmd += ["-with-report", f"{nextflow_out_dir}/run.report.html", "-with-timeline", f"{nextflow_out_dir}/run.timeline.html"]
   
    nf_cmd += ["-with-dag", f"{nextflow_out_dir}/run.dag.dot"]
    nf_cmd += ["-with-trace", f"{nextflow_out_dir}/run.trace.txt"]
    params_file = Path(nextflow_out_dir) / "run_params.yaml"
    nf_cmd += ["-params-file", str(params_file)]

    # Write params file
    if output:
        with open(params_file, 'w') as f:
            yaml.dump(task_params, f)
            f.flush()
    
    ##exit(1)

    if args.dry_run:
        print(" ".join(map(str, nf_cmd)))
    else:
        if args.verbosity >= VERBOSITY_VERBOSE:
            print(" ".join(map(str, nf_cmd)))
        resume_file = Path(nextflow_out_dir) / "resume.sh"
        with open(resume_file, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write(f"NXF_WORK={workdir} ")
            f.write(" ".join(map(str, nf_cmd)))
            f.write(" -resume")
            if os.environ.get('NXF_WORK'):
                f.write(" -work-dir " + os.environ['NXF_WORK'])
            f.write("\n")
        with open(Path(nextflow_out_dir) / "run_work_dir.txt", "wt") as f:
            f.write(workdir)
        try:
            env = os.environ.copy()
            env['NXF_WORK'] = workdir
            subprocess.run(nf_cmd, check=True, capture_output=(args.verbosity <= VERBOSITY_QUIET), text=True, env=env)
        except subprocess.CalledProcessError as e:
            print(e.stderr)
            print(f"To resume execution, run: sh {resume_file}")
            return 1
        collect_logs(nextflow_out_dir)

    if args.verbosity > VERBOSITY_QUIET:
        print(reads_msg, end='')
    if not args.dry_run and not args.stub_run:
        print_statistics(output)

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
