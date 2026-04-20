#!/usr/bin/env python
# requires pip install pyyaml
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
from urllib.parse import quote, urlparse
import json
import csv
import sqlite3
import glob
import hashlib
import tarfile
import urllib.request
import urllib.error
from typing import NamedTuple

import ssl
import urllib.parse
import uuid
import tempfile
from html.parser import HTMLParser

class SraMetadata(NamedTuple):
    """Metadata for an SRA run or synthetic equivalent for non-SRA reads."""
    run_accession: str
    sample_accession: str
    layout: str  # 'paired' or 'unpaired'
    read_count: str
    base_count: str
    avg_insert_size: str  # "NA" placeholder
    insert_size_stdev: str  # "NA" placeholder
    platform: str
    model: str
    experiment: str
    sra_study: str
    bio_sample: str
    project_id: str
    bio_project: str
    scientific_name: str
    tax_id: str
    release_date: str

def safe_urlopen(url, **kwargs):
    """Safe wrapper for urlopen that only allows http/https schemes."""
    parsed = urlparse(url)
    if parsed.scheme not in ('http', 'https'):
        raise ValueError(f"URL scheme '{parsed.scheme}' not allowed. Only http/https permitted.")
    return urlopen(url, **kwargs)


def safe_urlretrieve(url, filename):
    """Safe wrapper for urlretrieve that only allows http/https schemes."""
    parsed = urlparse(url)
    if parsed.scheme not in ('http', 'https'):
        raise ValueError(f"URL scheme '{parsed.scheme}' not allowed. Only http/https permitted.")
    return urllib.request.urlretrieve(url, filename)

# Requires pip install -r requirements.txt
import yaml

software_version = "0.5.1"
DEFAULT_DATA_VERSION = "current_1"

start_time = time.time()

VERBOSITY_DEFAULT=0
VERBOSITY_QUIET=-1
VERBOSITY_VERBOSE=1
verbosity = VERBOSITY_DEFAULT


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
ANIMALS = 33208   
## define gbdiv_inv as has ANIMALS but not VERTEBRATE
INSECTA = 50557
ARTHROPOD = 6656
VIRIDIPLANTAE = 33090
LEPIDOSAURS = 8504 
AMPHIBIANS = 8292 
ECHINODERMATA = 7586

# Reference organisms
H_SAPIENS = 9606
M_MUSCULUS = 10090
D_RERIO = 7955
D_MELANOGASTER = 7227
A_THALIANA = 3702
C_ELEGANS = 6239

BUSCO_DATA_URL = "https://busco-data.ezlab.org/v5/data"
FTP_EGAP_PROTOCOL = "https"
FTP_EGAP_SERVER = "ftp.ncbi.nlm.nih.gov"
FTP_EGAP_ROOT_PATH = "genomes/TOOLS/EGAP/support_data"
FTP_EGAP_ROOT = f"{FTP_EGAP_PROTOCOL}://{FTP_EGAP_SERVER}/{FTP_EGAP_ROOT_PATH}"
DATA_VERSION = os.environ.get("EGAPX_DATA_VERSION", DEFAULT_DATA_VERSION)
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
    parser.add_argument("-dv", "--data-version", help="Set support data version to specific value for run reproduction.", default="")
    parser.add_argument("-e", "--executor", help="Nextflow executor, one of docker, singularity, aws, or local (for NCBI internal use only). Uses corresponding Nextflow config file", default="local")
    parser.add_argument("-c", "--config-dir", help="Directory for executor config files, default is ./egapx_config. Can be also set as env EGAPX_CONFIG_DIR", default="")
    parser.add_argument("-w", "--workdir", help="Working directory for cloud executor", default="")
    parser.add_argument("--staging-dir", help="Staging subdirectory under workdir (for reproducible resume)", default="")
    parser.add_argument("-r", "--report", help="Report file prefix for report (.report.html) and timeline (.timeline.html) files, default is in output directory", default="")
    parser.add_argument("-n", "--dry-run", action="store_true", default=False)
    parser.add_argument("-st", "--stub-run", action="store_true", default=False)
    parser.add_argument("-so", "--summary-only", help="Print result statistics only if available, do not compute result", action="store_true", default=False)
    parser.add_argument("-lo", "--logs-only", help="Collect execution logs if available, put them in output directory, do not compute result", action="store_true", default=False)
    parser.add_argument("-ot", "--ortho-taxid", default=0, help="Taxid of reference data for orthology tasks")
    group = parser.add_argument_group('download')
    group.add_argument("-dl", "--download-only", help="Download external files to local storage, so that future runs can be isolated", action="store_true", default=False)
    group.add_argument("-dn", "--download-needed", dest='download_needed', help="Download only data needed by the current input into local cache", action="store_true", default=False)
    parser.add_argument("-lc", "--local-cache", help="Where to store the downloaded files", default="")
    parser.add_argument("-q", "--quiet", dest='verbosity', action='store_const', const=VERBOSITY_QUIET, default=VERBOSITY_DEFAULT)
    parser.add_argument("-v", "--verbose", dest='verbosity', action='store_const', const=VERBOSITY_VERBOSE, default=VERBOSITY_DEFAULT)
    parser.add_argument("-V", "--version", dest='report_version', help="Report software version", action='store_true', default=False)
    parser.add_argument("-fn", "--func_name", help="func_name", default="")
    parser.add_argument("--force", dest='force', help="Force execution despite of warnings", action='store_true', default=False)
    parser.add_argument("--ftp", dest='ftp', help="Enable FTP mode", action='store_true', default=False)
    
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


class HttpsDownloader:
    def __init__(self, verbosity=VERBOSITY_DEFAULT):
        self.base_url = f"https://{FTP_EGAP_SERVER}/{FTP_EGAP_ROOT_PATH}/"
        self.ssl_context = ssl.create_default_context()
        self.opener = urllib.request.build_opener(urllib.request.HTTPSHandler(context=self.ssl_context))
        self.host = f"https://{FTP_EGAP_SERVER}/"
        self.verbosity = verbosity

    def _url(self, path):
        """Build full URL from base_url + path (path should not start with /)."""
        if path.startswith("http://") or path.startswith("https://"):
            return path
        return f"{self.base_url}/{path.lstrip('/')}"
    
    def _head(self, url):
        """Return (headers dict, None) or (None, exception)."""
        req = urllib.request.Request(url, method="HEAD")
        try:
            with self.opener.open(req, timeout=30) as r:
                return (dict(r.headers), None)
        except Exception as e:
            return (None, e)

    def should_download_file(self, url, local_name, headers=None):
        """Return True if file should be downloaded (missing, or remote newer/larger)."""
        if headers is None:
            headers, err = self._head(url)
            if err or not headers:
                return True
            # can't HEAD, so try to download
        try:
            local_stat = os.stat(local_name)
        except FileNotFoundError:
            return True
        except NotADirectoryError:
            return True
        #local_mtime = datetime.datetime.fromtimestamp(local_stat.st_mtime)
        remote_size = headers.get("Content-Length")
        if remote_size is not None:
            remote_size = int(remote_size)
            if local_stat.st_size != remote_size:
                return True
        last_mod = headers.get("Last-Modified")
        if last_mod:
            try:
                from email.utils import parsedate_to_datetime
                remote_dt = parsedate_to_datetime(last_mod)

                if remote_dt.timestamp() > local_stat.st_mtime:
                    return True
            except Exception:
                pass
        return False

    def  download_file(self, url_or_path, local_path):
        """Download a single file from URL (or base_url + path) to local_path. Returns True on success."""
        url = self._url(url_or_path)
        os.makedirs(os.path.dirname(local_path) or ".", exist_ok=True)
        try:
            req = urllib.request.Request(url)
            with self.opener.open(req, timeout=60) as r:
                with open(local_path,"wb") as f:
                    f.write(r.read())
            return True
        except urllib.error.HTTPError as e:
            if e.code == 404:
                if self.verbosity >= VERBOSITY_DEFAULT:
                    print(f"FAILED 404: {url} -> {local_path}")
            else:
                if self.verbosity >= VERBOSITY_DEFAULT:
                    print(f"FAILED HTTP {e.code}: {url} -> {local_path}")
        except urllib.error.URLError as e:
            if self.verbosity >= VERBOSITY_DEFAULT:
                print(f"FAILED URL: {url} -> {local_path} - {e.reason}")
        except (BrokenPipeError, ConnectionError) as e:
            if self.verbosity >= VERBOSITY_DEFAULT:
                print(f"FAILED connection: {url} -> {local_path}")
            time.sleep(1)
            if self.verbosity >= VERBOSITY_DEFAULT:
                print("retrying...")
            return self.download_file(url_or_path, local_path)
        except Exception as e:
            print(f"FAILED: {url} -> {local_path} - {e}")
        return False

    def _parse_index_links(self, html, base_url):
        """Parse HTML directory index (Apache/nginx style) and return list of (href, is_dir). Skip parent dir."""
        links = []
        base = base_url if base_url.endswith("/") else base_url + "/"
        class IndexParser(HTMLParser):
            def handle_starttag(self,tag, attrs):
                if tag != "a":
                    return
                for k, v in attrs:
                    if k == "href" and v:
                        v = v.strip()
                        if not v or v in (".", "..")  or v.startswith("#") or v.lower().startswith("mailto:"):
                            return
                        full = urllib.parse.urljoin(base,v)
                        is_dir = v.endswith("/") or full.endswith("/")
                        links.append((full, is_dir))
                        break
        p = IndexParser()
        try:
            p.feed(html)
        except Exception:
            pass
        return links

    def _fetch_dir_listing(self,url):
        """Fetch URL; if HTML, return (True, links). Else return (False, None)."""
        try:
            req = urllib.request.Request(url)
            req.add_header("User-Agent","HttpsDownloader/1.0")
            with  self.opener.open(req,timeout=30) as r: 
                ct = r.headers.get("Content-Type", "")
                data = r.read()
                if "text/html" in ct:
                    try:
                        html = data.decode("utf-8", errors="replace")
                        links = self._parse_index_links(html, url)
                        return (True, links)
                    except Exception:
                        pass
        except Exception:
            pass
        return (False, None)

    def _is_parent_url(self, parent_url, child_url):
        """
        Checks if parent_url is a parent of child_url.
        """
        # 1. Parse both URLs into components
        parent_parts = urlparse(parent_url)
        child_parts = urlparse(child_url)

        # 2. Check if scheme and network location (domain) are the same
        if parent_parts.scheme != child_parts.scheme or parent_parts.netloc != child_parts.netloc:
            return False

        # 3. Use pathlib to compare the paths
        try:
            # Convert the URL paths to Path objects. We use PurePosixPath because URLs use forward slashes
            parent_path = Path(parent_parts.path).resolve()
            child_path = Path(child_parts.path).resolve()

            # is_relative_to checks if child_path can be expressed as a relative path from parent_path
            return child_path.is_relative_to(parent_path)
        except Exception:
            # Handle cases where path resolution might fail (e.g., non-existent paths on local machine)
            # If paths cannot be resolved, we can fall back to a simpler string-based check
            return child_parts.path.startswith(parent_parts.path)

    def download_dir(self, url_or_path, local_dir, only_if_newer=True):
        """Recursively download a directory. Works when the server returns an HTML directory index (e.g. Apache 'Index of').
        url_or_path: directory URL or path under base_url (should end with / for a directory).
        local_dir: local directory to write files into (created as needed).
        only_if_newer: if True, skip files that are already up to date (HEAD + mtime/size).

        """
        url = self._url(url_or_path)
        if not url.endswith("/"):
             url += "/"
        if self.verbosity > VERBOSITY_DEFAULT:
            print(f"Downloading directory {url} to {local_dir} (only_if_newer={only_if_newer})")
        os.makedirs(local_dir,exist_ok=True)
        is_index, links = self._fetch_dir_listing(url)
        if self.verbosity > VERBOSITY_DEFAULT:
            print(f"Fetched directory listing for {url}: is_index={is_index}, {len(links) if links else 0} links")
            print(f"Links: {links}")
        if not is_index or not links:
            return
        seen = set()
        for link_url, is_dir in links:
            # Normalize and skip duplicates / parent
            norm = urllib.parse.urlparse(link_url).path.rstrip("/")
            name = os.path.basename(norm) or os.path.basename(urllib.parse.urlparse(link_url).path.rstrip("/"))
            if self.verbosity > VERBOSITY_DEFAULT:
                print(f"Processing link: {link_url} (is_dir={is_dir}) -> name: {name}")
            if  not name or name in seen or name == "index.html" or name == "index.htm":
                continue
            if self._is_parent_url(link_url, url):
                if self.verbosity > VERBOSITY_DEFAULT:
                    print(f"Skipping parent URL: {link_url}")
                continue
            seen.add(name)
            local_path = os.path.join(local_dir, name)
            if self.verbosity > VERBOSITY_DEFAULT:
                print(f"Link {link_url} is {'directory' if is_dir else 'file'}, local path: {local_path}")
            if is_dir:
                # Ensure trailing slash for directory
                sub_url = link_url  if link_url.endswith("/") else link_url + "/"
                if self.verbosity > VERBOSITY_DEFAULT:
                    print(f"Recursively downloading directory {sub_url} to {local_path}")
                self.download_dir(sub_url, local_path, only_if_newer=only_if_newer)
            else:
                if only_if_newer and os.path.isfile(local_path):
                    if not self.should_download_file(link_url, local_path):
                        continue
                self.download_file(link_url, local_path)


def download_data_manifest(data_version):
    done = False
    while not done:
        done = True # Be optimistic :-)
        manifest_url = f"{FTP_EGAP_ROOT}/{data_version}.mft"
        manifest = safe_urlopen(manifest_url)
        manifest_list = []
        for line in manifest:
            line = line.decode("utf-8").strip()
            if not line or line[0] == '#':
                continue
            if line.endswith(".mft"):
                # Redirect and repeat
                data_version = line[:-len(".mft")]
                done = False
                break
            manifest_list.append(line)
    return manifest_list, data_version


def read_data_manifest(cache_dir, data_version):
    done = False
    err_msg = ""
    while not done:
        done = True # Be optimistic :-)
        manifest_list = []
        manifest_path = f"{cache_dir}/{data_version}.mft"
        if os.path.exists(manifest_path):
            with open(manifest_path, 'rt') as f:
                for line in f:
                    line = line.strip()
                    if not line or line[0] == '#':
                        continue
                    if line.endswith(".mft"):
                        # Redirect and repeat
                        data_version = line[:-len(".mft")]
                        done = False
                        break
                    manifest_list.append(line)
        else:
            err_msg = f"Data manifest '{manifest_path}' can not be read, please re-download data using -dl -lc {cache_dir}"
    return manifest_list, data_version, err_msg


def write_data_manifest_to_cache(cache_dir, manifest, data_version, effective_data_version):
    if not cache_dir:
        return
    manifest_path = os.path.join(cache_dir, effective_data_version + ".mft")
    with open(manifest_path, 'wt') as f:
        for line in manifest:
            f.write(f"{line}\n")
    if data_version != effective_data_version:
        master_manifest_path = os.path.join(cache_dir, data_version + ".mft")
        with open(master_manifest_path, 'wt') as f:
            f.write(effective_data_version + ".mft")


def download_egapx_data(cache_dir, data_version, use_ftp=False):
    manifest, effective_data_version = download_data_manifest(data_version)
    if data_version != effective_data_version:
        print(f"For replicating this run use --data-version {effective_data_version}")
    for line in manifest:
        print(f"Downloading {line}")
        if use_ftp:
            ftpd = FtpDownloader()
            ftpd.connect(FTP_EGAP_SERVER)
            ftpd.download_ftp_dir(f"{FTP_EGAP_ROOT_PATH}/{line}", os.path.join(cache_dir, line))
        else:
            # Implement HTTP download logic here if needed
            downloader = HttpsDownloader(verbosity=verbosity)
            downloader.download_dir(line, os.path.join(cache_dir, line))
    write_data_manifest_to_cache(cache_dir, manifest, data_version, effective_data_version)
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


def _get_uri_scheme(value):
    if not isinstance(value, str):
        return ''
    mo = re.match(r'^([a-zA-Z][a-zA-Z0-9+.-]*)://', value)
    return mo.group(1).lower() if mo else ''


def is_foreign_path(value, executor='local'):
    """A path is foreign if it needs staging before task execution."""
    if not isinstance(value, str) or value == '':
        return False
    scheme = _get_uri_scheme(value)
    if scheme in {'http', 'https', 'ftp'}:
        return True
    if executor == 'aws':
        return scheme != 's3'
    return False


def _join_staging_path(staging_root, rel_path):
    scheme = _get_uri_scheme(staging_root)
    if scheme:
        return f"{staging_root.rstrip('/')}/{rel_path.replace(os.sep, '/')}"
    return os.path.join(staging_root, rel_path)


def _build_staging_root(executor, workdir, staging_dir):
    stage_name = staging_dir or f"egapx_stage_{uuid.uuid4().hex[:12]}"
    work_root = workdir or 'work'
    work_scheme = _get_uri_scheme(work_root)
    if not work_scheme:
        work_root = os.path.abspath(work_root)
    return _join_staging_path(work_root, stage_name)


def _make_staging_target(source, staging_root):
    source_path = urlparse(source).path if _get_uri_scheme(source) else source
    base_name = os.path.basename(source_path.rstrip('/')) or 'downloaded_input'
    digest = hashlib.sha1(source.encode('utf-8')).hexdigest()[:12]
    return _join_staging_path(staging_root, os.path.join(digest, base_name))


def _convert_nested_paths(value, key, strict, leaf_converter):
    """Apply leaf_converter to scalar values while recursively walking dict/list."""
    if isinstance(value, dict):
        return {
            k: _convert_nested_paths(v, key, strict, leaf_converter)
            for k, v in value.items()
        }
    if isinstance(value, list):
        return [
            _convert_nested_paths(v, key, strict, leaf_converter)
            for v in value
        ]
    return leaf_converter(value, key, strict)


def _apply_path_converters(run_inputs, input_converter, output_converter):
    """Apply converters for known path-like keys in input and output sections."""
    input_root = run_inputs['input']
    for key in input_root:
        if key in path_inputs:
            strict = not key.endswith('reads')
            input_root[key] = input_converter(input_root[key], key, strict)

    if 'output' in run_inputs:
        run_inputs['output'] = output_converter(run_inputs['output'], 'output', False)


def _dedupe_transfer_plan(plan):
    """Keep deterministic order while removing duplicate source/target entries."""
    deduped = []
    seen = set()
    for item in plan:
        key = (item['source'], item['target'])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(item)
    return deduped


def convert_value(value, key, strict, executor='local', staging_root='', staging_plan=None, allow_staging=True):
    "Convert paths to absolute paths when possible, look into objects and lists"
    def _leaf(v, key_name, strict_mode):
        if not isinstance(v, str) or v == '':
            return v

        scheme = _get_uri_scheme(v)

        # Reuse locally cached foreign inputs when available (from -dn/-dl cache).
        if scheme in {'http', 'https', 'ftp', 's3'}:
            cache_dir = get_cache_dir()
            if cache_dir:
                cached_path = _find_cached_foreign_path(v, cache_dir)
                if cached_path:
                    return os.path.abspath(cached_path)

        foreign = allow_staging and is_foreign_path(v, executor)
        if foreign:
            # For local files we still need to resolve/check source before staging.
            source = v
            if not scheme:
                if not os.path.exists(v):
                    if strict_mode:
                        raise OSError(ENOENT, f"File for parameter '{key_name}' doesn't exist", v)
                    return v
                source = os.path.abspath(v)
            staged_target = _make_staging_target(source, staging_root)
            if staging_plan is not None:
                staging_plan.append({
                    'source': source,
                    'target': staged_target,
                    'key': key_name,
                    'scheme': scheme or 'file',
                    'foreign': True,
                })
            return staged_target

        if scheme:
            # Non-foreign URIs (for example s3:// on AWS) are kept as is.
            return v
        if not os.path.exists(v):
            if strict_mode:
                raise OSError(ENOENT, f"File for parameter '{key_name}' doesn't exist", v)
            return v
        # convert to absolute path
        return os.path.abspath(v)

    return _convert_nested_paths(value, key, strict, _leaf)


path_inputs = {'genome', 'hmm', 'softmask', 'reads_metadata', 'short_reads_metadata', 'long_reads_metadata', 'organelles',
               'proteins', 'proteins_trusted', 'additional_proteins', 'reads', 'short_reads', 'long_reads',
               'rnaseq_alignments', 'protein_alignments', 'ortho', 'reference_sets', 'prot_denylist',
               'name_cleanup_rules_file', 'gnomon_filtering_scores_file', 'busco_lineage_download',
               'cmsearch'}
def convert_paths(run_inputs):
    "Convert paths to absolute paths where appropriate"
    return convert_paths_for_executor(run_inputs)


def convert_paths_for_executor(run_inputs, executor='local', workdir='', staging_dir=''):
    """Convert paths and replace foreign inputs with staging targets.

    Returns a list of transfer plans to materialize in staging.
    """
    staging_root = _build_staging_root(executor, workdir, staging_dir)
    staging_plan = []
    _apply_path_converters(
        run_inputs,
        input_converter=lambda v, k, s: convert_value(
            v, k, s,
            executor=executor,
            staging_root=staging_root,
            staging_plan=staging_plan,
            allow_staging=True,
        ),
        output_converter=lambda v, k, s: convert_value(
            v, k, s,
            executor=executor,
            staging_root=staging_root,
            staging_plan=None,
            allow_staging=False,
        ),
    )
    return _dedupe_transfer_plan(staging_plan)


def stage_inputs_bulk(staging_plan, verbosity=VERBOSITY_DEFAULT):
    """Materialize staged inputs from local/http/ftp sources."""
    if not staging_plan:
        return True, ''

    msgs = []
    http_downloader = HttpsDownloader(verbosity=verbosity)
    ftp_downloaders = {}
    temp_root = tempfile.mkdtemp(prefix='egapx_stage_upload_')

    def _download_to_local_tmp(src, src_scheme):
        base_name = os.path.basename(urlparse(src).path.rstrip('/')) if src_scheme else os.path.basename(src)
        if not base_name:
            base_name = f"staged_{uuid.uuid4().hex}"
        local_tmp = os.path.join(temp_root, f"{uuid.uuid4().hex}_{base_name}")
        os.makedirs(os.path.dirname(local_tmp), exist_ok=True)
        if src_scheme in {'http', 'https'}:
            if not http_downloader.download_file(src, local_tmp):
                return '', f"Failed to download {src} -> {local_tmp}"
            return local_tmp, ''
        if src_scheme == 'ftp':
            parsed = urlparse(src)
            host = parsed.netloc
            if host not in ftp_downloaders:
                ftp_downloader = FtpDownloader()
                ftp_downloader.connect(host)
                ftp_downloaders[host] = ftp_downloader
            ftp_path = parsed.path.lstrip('/')
            if not ftp_path:
                return '', f"Invalid FTP source path: {src}"
            ftp_res = ftp_downloaders[host].download_ftp_file(ftp_path, local_tmp)
            if ftp_res is not True:
                return '', f"Failed to download FTP file {src} -> {local_tmp}"
            return local_tmp, ''
        return '', f"Unsupported temporary-download source scheme '{src_scheme}' for {src}"

    try:
        for item in staging_plan:
            src = item['source']
            dst = item['target']
            src_scheme = _get_uri_scheme(src)
            dst_scheme = _get_uri_scheme(dst)

            if dst_scheme == 's3':
                upload_src = src
                if src_scheme in {'http', 'https', 'ftp'}:
                    upload_src, err = _download_to_local_tmp(src, src_scheme)
                    if err:
                        return False, err
                elif src_scheme == '':
                    if not os.path.exists(src):
                        return False, f"Local source for staging doesn't exist: {src}"
                elif src_scheme != 's3':
                    return False, f"Unsupported staging source scheme '{src_scheme}' for {src}"

                try:
                    subprocess.run(
                        ["aws", "s3", "cp", upload_src, dst],
                        check=True,
                        stdout=None if verbosity >= VERBOSITY_VERBOSE else subprocess.DEVNULL,
                        stderr=None if verbosity >= VERBOSITY_VERBOSE else subprocess.DEVNULL,
                    )
                except FileNotFoundError:
                    return False, "AWS CLI not found in PATH, required for staging uploads to s3://"
                except subprocess.CalledProcessError:
                    return False, f"Failed to upload staged input {upload_src} -> {dst} using aws s3 cp"
                continue

            if dst_scheme in {'gs', 'az'}:
                return False, f"Staging target '{dst}' is not supported by this runner"

            os.makedirs(os.path.dirname(dst) or '.', exist_ok=True)

            if src_scheme in {'http', 'https'}:
                if not http_downloader.download_file(src, dst):
                    return False, f"Failed to download {src} -> {dst}"
            elif src_scheme == 'ftp':
                parsed = urlparse(src)
                host = parsed.netloc
                if host not in ftp_downloaders:
                    ftp_downloader = FtpDownloader()
                    ftp_downloader.connect(host)
                    ftp_downloaders[host] = ftp_downloader
                ftp_path = parsed.path.lstrip('/')
                if not ftp_path:
                    return False, f"Invalid FTP source path: {src}"
                ftp_res = ftp_downloaders[host].download_ftp_file(ftp_path, dst)
                if ftp_res is not True:
                    return False, f"Failed to download FTP file {src} -> {dst}"
            elif src_scheme == 's3':
                try:
                    subprocess.run(
                        ["aws", "s3", "cp", src, dst],
                        check=True,
                        stdout=None if verbosity >= VERBOSITY_VERBOSE else subprocess.DEVNULL,
                        stderr=None if verbosity >= VERBOSITY_VERBOSE else subprocess.DEVNULL,
                    )
                except FileNotFoundError:
                    return False, "AWS CLI not found in PATH, required for staging downloads from s3://"
                except subprocess.CalledProcessError:
                    return False, f"Failed to download staged input {src} -> {dst} using aws s3 cp"
            elif src_scheme == '':
                if not os.path.exists(src):
                    return False, f"Local source for staging doesn't exist: {src}"
                shutil.copy2(src, dst)
            else:
                return False, f"Unsupported staging source scheme '{src_scheme}' for {src}"
    finally:
        shutil.rmtree(temp_root, ignore_errors=True)

    return True, '\n'.join(msgs)


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

def prepare_reads(run_inputs, force=False, **kw):
    """Reformat reads input to be in 'fromFilePairs' format expected by egapx, i.e. [sample_id, [read1, read2]]
    Generate reads metadata file with minimal information - paired/unpaired and valid for existing libraries"""
    for_download = kw.get('for_download', False)
    res, msg, sras, nosra = True, '', list(), list()
    res1, msg1, sra1, nosra1 = prepare_reads_by_type(run_inputs, 'reads', 'short_reads', **kw)
    res, msg, sras, nosra = res & res1, msg + msg1, sras + sra1, nosra + nosra1 
    res1, msg1, sra1, nosra1 = prepare_reads_by_type(run_inputs, 'short_reads', None, **kw)
    res, msg, sras, nosra = res & res1, msg + msg1, sras + sra1, nosra + nosra1
    res1, msg1, sra1, nosra1  = prepare_reads_by_type(run_inputs, 'long_reads', None, **kw)
    res, msg, sras, nosra = res & res1, msg + msg1, sras + sra1, nosra + nosra1 

    sras = list(set(sras))
    # print(f'sras 424:   {sras}') 

    if 'output' in run_inputs and run_inputs['output'] and os.path.exists(run_inputs['output']):    
        outfile_path = Path(run_inputs['output'], 'sra_metadata.dat')
        save_sra_metadata_file(outfile_path, sras, nosra)
        run_inputs['input']['short_reads_metadata'] = str(outfile_path)
        run_inputs['input']['long_reads_metadata'] = str(outfile_path)

    # if force:    
    #     return res, msg, sras, nosra
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
            run_filter = None
            if key + '_ids' in run_inputs['input']:
                reads = run_inputs['input'][key + '_ids']
                run_filter = reads
                query = "[Accession] OR ".join(reads) + "[Accession]"
                if for_download:
                    del run_inputs['input'][key + '_ids']
                    run_inputs['input'][key + '_query'] = query
                    run_inputs['input'][key + '_filter'] = run_filter
                # Should we limit the number of runs similarly to query result?
            elif key + '_query' in run_inputs['input']:
                query = run_inputs['input'][key + '_query']
                real_query = True
            if query:
                sra_runs = sra_query(query, run_filter)
                # print(f"sra_runs: {sra_runs}")
                sra_runs_ids = list()
                paired = 0
                unpaired = 0
                unidentified = 0
                for rec in sra_runs:
                    sra_runs_ids.append(rec.run_accession)
                    if rec.layout == 'paired':
                        paired += 1
                    elif rec.layout == 'unpaired':
                        unpaired += 1
                    else:
                        unidentified += 1
                if unidentified > 0:
                    msg += 'ERROR: incorrect metadata format returned'
                    res = False
                else:
                    msg += f"Detected {unpaired} single-end runs, {paired} paired-end runs for query '{key}'\n"
                if real_query and not for_download:
                    del run_inputs['input'][key + '_query']
                    run_inputs['input'][key + '_ids'] = sra_runs_ids
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
                
    return res, msg


def is_sra(run_name):
    return re.fullmatch(r'([DES]RR[0-9]+)', run_name) is not None


def _pop_reads_input(run_inputs, reads_type):
    if reads_type in run_inputs['input']:
        reads = run_inputs['input'][reads_type]
        del run_inputs['input'][reads_type]
        return True, reads
        
    reads_type_ids = f"{reads_type}_ids"
    if reads_type_ids in run_inputs['input']:
        reads = run_inputs['input'][reads_type_ids]
        del run_inputs['input'][reads_type_ids]
        return True, reads
        
    reads_type_query = f"{reads_type}_query"
    if reads_type_query in run_inputs['input']:
        reads = run_inputs['input'][reads_type_query]
        del run_inputs['input'][reads_type_query]
        return True, reads
        
    return False, None


def _parse_reads_table_file(filename):
    prefixes = defaultdict(list)
    sra_to_return = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            mo = re.match(r'([^ \t]+)[ \t]+(.*)', line)
            if not mo:
                continue
            sample = mo.group(1)
            files_part = mo.group(2).strip()
            files = files_part.split()
            prefixes[sample].extend(files)
    return prefixes, sra_to_return


def _handle_query_string(reads_query, reads_type_write, run_inputs, run_filter=None):
    metadata_map, file_map = sra_query_cached(reads_query, run_filter)
    if file_map:
        return True, defaultdict(list, file_map), list(metadata_map.keys()), metadata_map, ''
    run_inputs['input'][f'{reads_type_write}_query'] = reads_query
    return False, defaultdict(list), list(metadata_map.keys()), metadata_map, ''


def _collect_from_flat_list_entry(rf, prefixes, sra_to_return, msg):
    if isinstance(rf, str):
        name = Path(rf).name
        mo = re.match(r'([^._]+)', name)
        if mo and mo.group(1) != name:
            run_name = mo.group(1)
            if is_sra(run_name):
                sra_to_return.append(run_name)
            prefixes[run_name].append(rf)
            return True, msg
        elif is_sra(rf):
            sra_to_return.append(rf)
            fnames = get_cached_sra_files(rf)
            if fnames:
                prefixes[rf].extend(fnames)
            return True, msg
        else:
            msg += f"Invalid read input {rf} - neither a file name nor SRA run name\n"
            return False, msg
    msg += f"Invalid read input type for {rf}\n"
    return False, msg


def _infer_sample_id_from_file_list(file_list):
    # Ensure all elements are strings
    if not all(isinstance(x, str) for x in file_list):
        raise ValueError("All entries in file_list must be strings")
    names = [re.match(r'([^.]+)', Path(x).name).group(1) for x in file_list]
    names.sort()
    if len(names) == 1:
        sample_id = names[0]
    else:
        limit = min(len(names[0]), len(names[-1]))
        i = 0
        while i < limit and names[0][i] == names[-1][i]:
            i += 1
        sample_id = names[0][0:i]
    while sample_id and sample_id[-1] in '._':
        sample_id = sample_id[:-1]
    return sample_id

def _process_list_reads(reads: list, prefixes: defaultdict(list), sra_to_return: list) -> (bool, str):
    """ Process reads input given as a list
    :param reads: list of reads entries
    :param prefixes: map of run_name to list of read files
    :param sra_to_return: list of SRA run names to return
    :return: success flag, error message
    """
    msg = ''
    for rf in reads:
        if isinstance(rf, str):
            ok, msg = _collect_from_flat_list_entry(rf, prefixes, sra_to_return, msg)
            if not ok:
                return False, msg
        elif isinstance(rf, list):
            if len(rf) == 2 and isinstance(rf[0], str) and isinstance(rf[1], list):
                run_name, run_files = rf
                if run_name in prefixes:
                    msg += f"Run name {run_name} used multiple times\n"
                    return False, msg
                for r in run_files:
                    if not isinstance(r, str):
                        msg += f"Invalid read input {r}\n"
                        return False, msg
                    prefixes[run_name].append(r)
                if is_sra(run_name):
                    sra_to_return.append(run_name)
            else:
                try:
                    sample_id = _infer_sample_id_from_file_list(rf)
                except Exception as e:
                    msg += f"Invalid read input in file list: {rf} ({e})\n"
                    return False, msg
                prefixes[sample_id].extend(rf)
        else:
            msg += f"Unsupported reads entry type: {type(rf)}\n"
            return False, msg
    return True, msg


def _generate_fake_metadata(prefixes: defaultdict[str, list[str]], sra_metadata_map: dict[str, SraMetadata]) -> tuple[list[SraMetadata], str]:
    """ Generate minimal SRA metadata for non-SRA read inputs
    :param prefixes: map of run_name to list of read files
    :param sra_metadata_map: map of SRA run_name to existing metadata records
    :return: list of generated metadata records, error message
    """
    not_sra = []
    error_msg = ''
    for run_name, files in prefixes.items():
        mos = [re.search(r'([DES]RR[0-9]+)', x) for x in files]
        srrs = list({mo.group(1) for mo in mos if mo})
        if len(srrs) > 1:
            error_msg += f"Error: Run set {run_name} has multiple run ids in read files: {srrs}\n"
            continue
        run_name_resolved = srrs[0] if srrs else run_name
        if is_sra(run_name_resolved):
            paired = 'paired' if len(files) == 2 else 'unpaired'
            if run_name_resolved in sra_metadata_map:
                out_rec = list(sra_metadata_map[run_name_resolved])
                if len(out_rec) == 13:
                    out_rec.append('NA') # Insert 'NA' for bio_project
                    out_rec.append('NA') # Insert 'NA' for scientific_name
                    out_rec.append('0')  # Insert '0' for tax_id
                    out_rec.append('NA') # Insert 'NA' for release_date
            else:
                out_rec = [run_name_resolved, 'NA', paired, '2', '2', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA',
                           'SAMN_' + run_name_resolved, '0', 'NA', 'NA', '0', 'NA']
            not_sra.append(SraMetadata(*out_rec))
    return not_sra, error_msg


def _validate_file_counts(prefixes):
    for run_name, files in prefixes.items():
        if len(files) > 2:
            return False, f"More than 2 files for run {run_name}\n"
    return True, ''


def prepare_reads_by_type(run_inputs, reads_type, reads_type_write=None, **kw):
    for_download = kw.get('for_download', False)
    fake_metadata = kw.get('fake_metadata', not for_download)
    # print(f"Preparing reads for type '{reads_type}' with fake_metadata={fake_metadata}")
    res = True
    msg = ''
    sra_metadata_map = {}

    present, reads = _pop_reads_input(run_inputs, reads_type)
    # print(f"Reads input for type '{reads_type}': {reads}")
    if not present:
        return res, msg, [], []
    if reads_type_write is None:
        reads_type_write = reads_type

    prefixes = defaultdict(list)
    sra_to_return = []
    has_files = False

    if isinstance(reads, str):
        if os.path.exists(reads):
            has_files = True
            prefixes, sra_to_return = _parse_reads_table_file(reads)
        else:
            has_files, prefixes, sra_to_return, sra_metadata_map, add_msg = _handle_query_string(
                reads, reads_type_write, run_inputs
            )
            msg += add_msg
            if not has_files:
                return res, msg, sra_to_return, []
    elif isinstance(reads, list):
        ok, list_msg = _process_list_reads(reads, prefixes, sra_to_return)
        msg += list_msg
        if not ok:
            return False, msg, [], []
        has_files = any(prefixes.values())
    else:
        msg += f"Unsupported reads input type: {type(reads)}\n"
        return False, msg, [], []

    not_sra = []
    if has_files:
        if fake_metadata:
            not_sra, gen_err = _generate_fake_metadata(prefixes, sra_metadata_map)
            if gen_err:
                msg += gen_err
                return False, msg, sra_to_return, not_sra
        ok, vmsg = _validate_file_counts(prefixes)
        msg += vmsg
        if not ok:
            return False, msg, sra_to_return, not_sra
        run_inputs['input'][reads_type_write] = [[k, v] for k, v in prefixes.items()]
    elif reads:
        run_inputs['input'][f'{reads_type_write}_ids'] = reads

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
    
    if 'lineage_taxids' not in inputs:
        lineage = get_lineage(taxid)
        inputs['lineage_taxids'] = ','.join(map(str,lineage))

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
        aligner_name = inputs.get('protein_aligner_name', 'miniprot')
        # Use best N taxons from the best bin for selecting proteins from the set,
        # depends on the aligner algorithm user and whether it is set explicitly
        if aligner_name == "prosplign":
            best_n = 5
        else:
            best_n = 10
        best_n = int(inputs.get('proteins_best_n_orgs', best_n))
        
        proteins, trusted, taxid_filter = get_closest_protein_bag(
            taxid, best_n, set(map(int, inputs.get('proteins_deny_taxids', []))))
        if not proteins:
            # Proteins are not specified explicitly and not found by taxid
            print(f"ERROR: Proteins are not found for tax id {inputs['taxid']}")
            print("  This organism is not supported by current NCBI protein database")
            print("  You can specify proteins manually using 'proteins' parameter")
            return False
        inputs['proteins'] = proteins
        if trusted:
            inputs['proteins_trusted'] = trusted
        if taxid_filter:
            inputs['proteins_filter_taxons'] = taxid_filter

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

    if inputs.get('cmsearch', {}).get('enabled'):
        inputs['cmsearch'] = {'files': [get_file_path('cmsearch', f) for f in "Rfam.seed rfam1410.cm rfam1410_amendments.xml".split()] }
        inputs['cmsearch']['enabled'] = True
        if len(inputs['cmsearch']['files']) != 3:
            print(f"ERROR: Expected 3 cmsearch files (Rfam.seed, rfam1410.cm, rfam1410_amendments.xml), but found {len(inputs['cmsearch']['files'])}.")
            print("  Check that the cmsearch data files are properly configured.")
            return False
    else:
        inputs['cmsearch'] = {'enabled':False}

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
def load_version_map(data_version):
    """
    Load version map from manifest file either from local data directory or from NCBI and cache it.
    """
    global data_version_cache
    cache_dir = get_cache_dir()
    if cache_dir:
        manifest, effective_data_version, err_msg = read_data_manifest(cache_dir, data_version)
        if err_msg:
            return effective_data_version, err_msg
    else:
        manifest, effective_data_version = download_data_manifest(data_version)
    for line in manifest:
        parts = line.split('/')
        if len(parts) == 2:
            data_version_cache[parts[0]] = parts[1]
    return effective_data_version, ''


def get_version_map():
    return data_version_cache


def get_versioned_path(subsystem, filename) -> str:
    """
    Retrieve the versioned path for a given subsystem and filename.

    Returns the path including the version if found; otherwise, returns the default path.
    """
    data_version_map = get_version_map()
    if subsystem not in data_version_map:
        return os.path.join(subsystem, filename)
    version = data_version_map[subsystem]
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
        # Create directory and copy user config files there
        from_dir = Path(script_directory) / 'assets' / 'config' / 'user'
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
        #mo = re.search(r"process.+(memory|cpus) *=", config_txt)
        #if not mo:
        #    default_configs.append("process_resources.config")
    
    # Add mandatory configs
    config_files = [str(config_file)]
    for cf in default_configs:
        config_files.insert(0, os.path.join(script_directory, "assets/config", cf))
    return ",".join(config_files)


lineage_cache = {}
taxname_cache = {}

def resilient_urlopen(url, max_retries=3, retry_delay=1.0):
    """
    Wrapper for safe_urlopen with retry logic for handling intermittent failures.
    
    Args:
        url: URL to open
        max_retries: Maximum number of retry attempts (default: 3)
        retry_delay: Initial delay between retries in seconds (default: 1.0)
                     Uses exponential backoff
    
    Returns:
        Response object from urlopen
    
    Raises:
        URLError or HTTPError if all retries fail or on non-retryable errors
    """
    delay = retry_delay
    for attempt in range(max_retries):
        try:
            return safe_urlopen(url)
        except (urllib.error.URLError, urllib.error.HTTPError) as e:
            # Check if it's a retryable error
            is_retryable = False
            if isinstance(e, urllib.error.HTTPError):
                # Retry on 503 (service unavailable) and 429 (rate limited)
                is_retryable = e.code in (503, 429)
            else:
                # Retry on network errors (URLError)
                is_retryable = True
            
            if is_retryable and attempt < max_retries - 1:
                # Not the last attempt, wait and retry
                time.sleep(delay)
                delay *= 2  # Exponential backoff
            else:
                # Last attempt or non-retryable error, propagate
                raise

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
    taxon_json_file = resilient_urlopen(dataset_taxonomy_url+str(taxid))
    taxon = json.load(taxon_json_file)["taxonomy_nodes"][0]
    if 'errors' in taxon:
        lineage_cache[taxid] = [], {}
        return [], {}
    # print(f"taxid: {taxid}, taxon: {taxon}")
    lineage = taxon["taxonomy"].get("lineage", [])
    lineage.append(taxon["taxonomy"]["tax_id"])

    ranks_file = resilient_urlopen(dataset_taxonomy_url+",".join(map(str, lineage)))
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


def get_tax_name(taxid) -> str:
    """
    Returns the taxonomic name for a given taxid.
    NB: the name cache is valid only AFTER getting the lineage for the taxid or one of its descendants
    """
    global taxname_cache
    return taxname_cache.get(taxid, '')


def get_tax_file(subsystem, tax_path):
    vfn = get_versioned_path(subsystem, tax_path)
    cache_dir = get_cache_dir()
    taxids_url = f"{FTP_EGAP_ROOT}/{vfn}"
    taxids_path = os.path.join(cache_dir, vfn) if cache_dir else ""

    if cache_dir and os.path.exists(taxids_path):
        with open(taxids_path, "rb") as r:
            return r.readlines()

    try:
        response = safe_urlopen(taxids_url)
        contents = response.read()
        if cache_dir:
            os.makedirs(os.path.dirname(taxids_path), exist_ok=True)
            with open(taxids_path, "wb") as w:
                w.write(contents)
        return contents.splitlines(keepends=True)
    except urllib.error.HTTPError:
        raise FileNotFoundError


def read_versioned_file(subsystem, fname):
    vfn = get_versioned_path(subsystem, fname)
    cache_dir = get_cache_dir()
    taxids_url = f"{FTP_EGAP_ROOT}/{vfn}"
    file_path = os.path.join(cache_dir, vfn) if cache_dir else ""

    if cache_dir and os.path.exists(file_path):
        return Path(file_path).read_text(encoding='utf-8', errors='ignore')

    try:
        response = safe_urlopen(taxids_url)
        raw = response.read()
        encoding = response.headers.get_content_charset() or 'utf-8'
        contents = raw.decode(encoding)
        if cache_dir:
            os.makedirs(os.path.dirname(file_path), exist_ok=True)
            with open(file_path, 'wb') as w:
                w.write(raw)
        return contents
    except urllib.error.HTTPError:
        raise FileNotFoundError


def check_supported_taxid(taxid):
    if not taxid:
        return False, "No taxid provided"

    lineage = get_lineage(taxid)
    supported_taxa = yaml.safe_load(read_versioned_file("misc", "supported_taxa.yaml"))
    names = []
    reqs = set()
    for name, taxid in supported_taxa:
        names.append(f"{name}({taxid})")
        reqs.add(taxid)
    lineage_check = any(l in lineage for l in reqs)

    #print(f'input tax: {taxid} , lineage: {lineage} , contains: { lineage_check } ')
    if lineage_check:
        return True, ''
    supported_names = ", ".join(names[:-1]) + ", or " + names[-1]
    return False, f"ERROR: Input taxid not within supported range. Must be under {supported_names} according to NCBI Taxonomy (https://www.ncbi.nlm.nih.gov/taxonomy)"


def compare_lineages(ln1, ln2):
    min_common_len = min(len(ln1), len(ln2))
    for i in range(min_common_len):
        if ln1[i] != ln2[i]:
            return i - 1
    return min_common_len - 1


def get_protein_filter(
    bin_taxid : int,     # taxid of selected bin to filter from 
    lineage : List[int], # root-first taxid list for main organism
    best_n : int,        # number of organisms to select from bin, mandatory not counted
    deny_set : set[int]  # set of taxids not to be included
) -> List[str]:
    # print(f"get_protein_filter({bin_taxid}, {lineage}, {best_n}, {deny_set})")
    # Find best_n organisms in the selected bin
    org_matches = []
    # To preserve original order in file for organisms with the same score
    ordinal = 0
    max_pos = len(lineage)
    bin_file = get_tax_file("target_proteins", f"{bin_taxid}.assembly.taxid.list")
    taxid_name_map = {}
    layed_out_records = []
    for line in bin_file:
        line = line.decode('utf-8').strip()
        if not line or line[0] == '#':
            continue
        parts = line.split('\t')
        # print(parts)
        # taxon, name, lineage (semicolon-separated)
        taxon, name, lineage_str = parts
        t = int(taxon)
        if t in deny_set:
            continue
        org_lineage = list(map(int, lineage_str.split(';')))
        org_lineage.append(t)
        score = max_pos - compare_lineages(lineage, org_lineage)
        if t == 10090:
            layed_out_records.append((score, ordinal, t, name))
        else:
            org_matches.append((score, ordinal, t, name))
        taxid_name_map[t] = name
        ordinal += 1
    mandatory_taxids = {H_SAPIENS, M_MUSCULUS, D_RERIO, D_MELANOGASTER, A_THALIANA, C_ELEGANS}
    # To make sure that both H_SAPIENS and M_MUSCULUS are not in the list
    # we add M_MUSCULUS only if we did not see H_SAPIENS
    if H_SAPIENS not in taxid_name_map:
        org_matches += layed_out_records
    else:
        # Exclude M_MUSCULUS if H_SAPIENS is present
        mandatory_taxids -= {M_MUSCULUS}
    org_matches.sort()
    # print(f"org_matches: {org_matches}")
    # Select best_n first entries from the sorted org_matches,
    # organisms in deny_set were discarded when reading the file.
    protein_filter_set = { rec[2] for rec in org_matches[:best_n] }
    # Extend with organisms from org_matches and mandatory_taxids,
    # excluding 10090 if 9606 is already present
    mandatory_taxids &= taxid_name_map.keys()
    protein_filter_set |= mandatory_taxids

    return [taxid_name_map[t] for t in protein_filter_set]


def get_closest_protein_bag(taxid, best_n, deny_set={}):
    """
    Using taxid as a source find the closest protein bag.
    In the protein bag find best_n organisms but exclude organisms in deny_set.

    """
    # print(f"get_closest_protein_bag({taxid}, {best_n}, {deny_set})")
    no_result = '', '', []
    if not taxid:
        return no_result

    new_format_protein_bins = False
    taxid_list = []
    taxids_file = get_tax_file("target_proteins", "taxid.list")
    for line in taxids_file:
        line = line.decode("utf-8").strip()
        if not line or line[0] == '#':
            continue
        parts = line.split('\t')
        t = parts[0]
        taxid_list.append(int(t))
    # Check for new format
    try:
        # Try first file in main taxid.list
        get_tax_file("target_proteins", f"{taxid_list[0]}.assembly.taxid.list")
        new_format_protein_bins = True
    except FileNotFoundError:
        # Old format
        new_format_protein_bins = False

    lineage = get_lineage(taxid)

    best_taxid = 0
    best_score = 0
    for t in taxid_list:
        try:
            pos = lineage.index(t)
        except ValueError:
            continue
        if pos > best_score:
            best_score = pos
            best_taxid = t

    if best_score == 0:
        return no_result
    if new_format_protein_bins and best_n > 0:
        protein_filter = get_protein_filter(best_taxid, lineage, best_n, deny_set)
    else:
        protein_filter = []

    # print(f"protein_filter: {protein_filter}, protein_filter_set: {protein_filter_set}")
    if new_format_protein_bins:
        protein_file = get_file_path("target_proteins", f"{best_taxid}.proteins.faa.gz")
    else:
        protein_file = get_file_path("target_proteins", f"{best_taxid}.faa.gz")
    trusted_file = get_file_path("target_proteins", f"{best_taxid}.trusted_proteins.gi")
    return protein_file, trusted_file, protein_filter


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
    best_lineage = []
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
    # best_lineage = None
    # best_good_match = False
    best_taxid = None
    best_score = 0
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
            # best_lineage = l
            # best_good_match = last_good_match

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
        safe_urlretrieve(busco_list_url, busco_list_file)


def check_busco_lineage(cache_dir, lineage):
    busco_downloads_dir = os.path.join(cache_dir, "busco_downloads")
    lineages_path = os.path.join(busco_downloads_dir, "lineages")
    return check_busco_lineage_in_lineages_path(lineages_path, lineage)


def check_busco_lineage_in_lineages_path(lineages_path, lineage):
    lineage_path = os.path.join(lineages_path, lineage)
    if not os.path.isdir(lineage_path) or not os.listdir(lineage_path):
        return False

    # Some externally prepared BUSCO bundles contain only lineages/<lineage>
    # without a companion <lineage>.list file. Accept these as valid.
    lineage_list = os.path.join(lineages_path, f"{lineage}.list")
    if os.path.exists(lineage_list):
        for line in open(lineage_list, "rt"):
            line = line.strip()
            if not line:
                continue
            if not os.path.exists(os.path.join(lineages_path, line)):
                return False
    return True


def resolve_busco_lineage_download_path(busco_input_path, lineage):
    """Resolve and validate BUSCO lineage input path.

    Supported forms:
    - BUSCO downloads root containing file_versions.tsv and lineages/
    - lineages/ directory containing <lineage>/
    - direct path to the required lineage directory
    """
    if not isinstance(busco_input_path, str) or not busco_input_path.strip():
        return False, '', "Parameter 'busco_lineage_download' must be a non-empty path"

    candidate = os.path.abspath(busco_input_path)
    if not os.path.isdir(candidate):
        return False, '', f"BUSCO path does not exist or is not a directory: {candidate}"

    # Case 1: root directory with file_versions.tsv and lineages/
    root_list = os.path.join(candidate, "file_versions.tsv")
    root_lineages = os.path.join(candidate, "lineages")
    if os.path.isfile(root_list) and os.path.isdir(root_lineages):
        if check_busco_lineage_in_lineages_path(root_lineages, lineage):
            return True, os.path.join(root_lineages, lineage), ''
        return False, '', f"BUSCO lineage '{lineage}' not found or invalid under {root_lineages}"

    # Case 2: lineages directory
    if check_busco_lineage_in_lineages_path(candidate, lineage):
        return True, os.path.join(candidate, lineage), ''

    # Case 3: direct lineage directory
    if os.path.basename(candidate.rstrip('/')) == lineage and os.listdir(candidate):
        return True, candidate, ''

    return False, '', (
        f"BUSCO path {candidate} is neither a BUSCO downloads root, lineages directory, "
        f"nor a valid lineage directory for '{lineage}'"
    )


def download_busco_lineage(cache_dir, lineage):
    if check_busco_lineage(cache_dir, lineage):
        return
    busco_downloads_dir = os.path.join(cache_dir, "busco_downloads")
    lineages_path = os.path.join(busco_downloads_dir, "lineages")
#    if os.path.isdir(lineage_path) and os.listdir(lineage_path):
#        return
    lineage_set = get_busco_lineage_set()
    lineage_date = lineage_set.get(lineage)
    if not lineage_date:
        raise Exception(f"Lineage {lineage} not found")
    # Download the lineage tar.gz and unpack it
    busco_lineage_url = f"{BUSCO_DATA_URL}/lineages/{lineage}.{lineage_date}.tar.gz"
    busco_lineage_file = os.path.join(busco_downloads_dir, f"{lineage}.tar.gz")
    os.makedirs(lineages_path, exist_ok=True)
    safe_urlretrieve(busco_lineage_url, busco_lineage_file)
    with tarfile.open(busco_lineage_file) as tar:
        with open(os.path.join(lineages_path, f"{lineage}.list"), "wt") as tar_list:
            tar_list.write("\n".join(tar.getnames()))
        tar.extractall(path=lineages_path)
    os.remove(busco_lineage_file)


def ensure_busco_lineage_for_inputs(cache_dir, inputs):
    """Ensure BUSCO lineage data is available for run inputs unless pre-specified."""
    # If busco_lineage_download is set, it means lineage data was provided by user.
    if inputs.get('busco_lineage_download'):
        return

    download_busco_lineage_set(cache_dir)
    busco_lineage = inputs.get('busco_lineage') or get_busco_lineage(inputs['taxid'])
    print(f"Downloading BUSCO lineage {busco_lineage}")
    download_busco_lineage(cache_dir, busco_lineage)


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
        return { x[0]: x[1] for x in map(lambda x: x.decode("utf-8").split('\t'), safe_urlopen(busco_list_url).readlines())}


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
            _, task_hash, _, name = parts[:4]
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
            hash_parts = task_hash.split("/")
            if workdir.startswith("s3://"):
                # AWS S3
                task_workdir_stub = os.path.join(workdir, task_hash)
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
                lines = os.listdir(task_workdir_stub)
                # print(lines)
                for line in lines:
                    if line.startswith(hash_parts[1]):
                        task_workdir = os.path.join(task_workdir_stub, line)
                        break

                for src_file in glob.glob(os.path.join(task_workdir, ".command.*")):
                    shutil.copy(src_file, task_dst)

            logs = glob.glob(os.path.join(task_dst, ".command.*"))
            # Make logs visible
            for log in logs:
                os.rename(log, os.path.join(task_dst, os.path.basename(log)[1:]))
            print_log_values(outdir, task_key, g_report_effective_values)
            # print(f"collected from {task_hash}, task {name} {task_workdir} to {task_dst}")
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
    # print(f"sra query {sra_potential_query}")
    time.sleep(1)
    biomol_str = " AND biomol_transcript[properties] "
     
    if 'biomol_transcript' not in sra_potential_query:
        sra_potential_query += biomol_str
   
    esearch = safe_urlopen(esearch_url+quote(sra_potential_query.encode('utf-8')))
    esearch_json = json.load(esearch)
    uids = esearch_json["esearchresult"]["idlist"]
    return uids
    ##uids = ','.join(esearch_json["esearchresult"]["idlist"])
    ##print(esearch_json)


def sra_metadata_query(sra_uids_list: list[str], run_filter: set[str] | None = None) -> list[SraMetadata]:
    """Query NCBI for SRA run metadata.
    
    Args:
        sra_uids_list: List of SRA UIDs to query
        run_filter: Optional set of run accessions to filter results
        
    Returns:
        List of SraMetadata records
    """
    uids = ','.join(sra_uids_list)
    runinfo = safe_urlopen(runinfo_url+uids)
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
    ##            "Platform type", "Model", "SRA Experiment accession", "SRA Study accession", "Biosample accession", "Bioproject ID",
    ##            "Bioproject Accession", "Scientific name", "TaxID", "Release date"]
    
    records: dict[str, SraMetadata] = {}
    for parts in csv.reader(body, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True):
        if len(parts) == 0:
            continue
        if "Run" not in keypos:
            continue
        run = parts[keypos["Run"]]
        if run_filter is not None and run not in run_filter:
            continue
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
        records[run] = SraMetadata(*record)
        # print(TAB.join(record))
    return records.values()


g_sra_to_file_map: dict[str, list[str]] = {}
g_query_to_accessions_map: dict[str, list[str]] = {}
g_sra_metadata_map: dict[str, SraMetadata] = {}

def read_cache() -> tuple[dict[str, list[str]], dict[str, list[str]], dict[str, SraMetadata]]:
    """Read cached SRA data from disk.
    
    Returns:
        tuple of:
            - sra_to_file_map: run accession -> list of file paths
            - query_to_accessions_map: query hash -> list of run accessions
            - sra_metadata_map: run accession -> SraMetadata
    """
    global g_sra_to_file_map
    global g_query_to_accessions_map
    global g_sra_metadata_map
    cache_dir = get_cache_dir()
    if cache_dir:
        if not g_sra_to_file_map:
            sra_runs_file = os.path.join(cache_dir, SRA_DOWNLOAD_FOLDER, SRA_RUNS_FILE)
            if os.path.isfile(sra_runs_file):
                g_sra_to_file_map = yaml.safe_load(open(sra_runs_file, 'r'))
        if not g_query_to_accessions_map:
            sra_queries_file = os.path.join(cache_dir, SRA_DOWNLOAD_FOLDER, SRA_QUERIES_FILE)
            if os.path.isfile(sra_queries_file):
                g_query_to_accessions_map = yaml.safe_load(open(sra_queries_file, 'r'))
        if not g_sra_metadata_map:
                sra_metadata_file = os.path.join(cache_dir, SRA_DOWNLOAD_FOLDER, SRA_METADATA_FILE)
                if os.path.isfile(sra_metadata_file):
                    raw_metadata = yaml.safe_load(open(sra_metadata_file, 'r'))
                    if raw_metadata:
                        for k, v in raw_metadata.items():
                            if isinstance(v, list) and len(v) == 17:
                                g_sra_metadata_map[k] = SraMetadata(*v)
                            elif isinstance(v, dict):
                                g_sra_metadata_map[k] = SraMetadata(**v)
    return g_sra_to_file_map, g_query_to_accessions_map, g_sra_metadata_map


def sra_query_cached(query, run_filter=None) -> tuple[dict[str, SraMetadata], dict[str, list[str]]]:
    """Query SRA with caching support.
    
    Args:
        query: Entrez query string
        run_filter: Optional list of run accessions to filter results
        
    Returns:
        tuple of:
            - metadata_map: run accession -> SraMetadata
            - file_map: run accession -> list of cached file paths
    """
    sra_to_file_map, query_to_accessions_map, sra_metadata_map = read_cache()
    metadata_map: dict[str, SraMetadata] = {}
    file_map: dict[str, list[str]] = {}

    query_hash = hashlib.sha1(query.encode()).hexdigest()
    if query_to_accessions_map and query_hash in query_to_accessions_map:
        accessions = query_to_accessions_map.get(query_hash, {})
        if accessions:
            runs = accessions.get('runs', [])
            # TODO: should we apply filter logic here as well?
            for run in runs:
                fnames = sra_to_file_map.get(run)
                if fnames:
                    file_map[run] = fnames
                if run in sra_metadata_map:
                    metadata_map[run] = sra_metadata_map[run]
    else:
        metadata = sra_query(query, run_filter)
        if metadata:
            for rec in metadata:
                metadata_map[rec.run_accession] = rec
    
    return metadata_map, file_map


def get_cached_sra_files(run_name):
    sra_to_file_map, _, _ = read_cache()
    return sra_to_file_map.get(run_name)


def get_sra_metadata(runs) -> dict[str, SraMetadata]:
    """Get metadata for a list of SRA runs.
    
    Args:
        runs: List of SRA run accessions
        
    Returns:
        Map of run accession -> SraMetadata
    """
    metadata_map: dict[str, SraMetadata] = {}
    _, _, sra_metadata_map = read_cache()
    if sra_metadata_map:
        for run in runs:
            if run in sra_metadata_map:
                metadata_map[run] = sra_metadata_map[run]
    else:
        metadata = sra_query(" OR ".join(runs), runs)
        if metadata:
            for rec in metadata:
                metadata_map[rec.run_accession] = rec

    return metadata_map


def sra_query(query, run_filter=None) -> list[SraMetadata]:
    """Query SRA and return metadata for matching runs.
    
    Args:
        query: Entrez query string
        run_filter: Optional list of run accessions to filter results
        
    Returns:
        List of SraMetadata records
    """
    global verbosity
    if verbosity >= VERBOSITY_VERBOSE:
        print(f"sra_query: query '{query}', run_filter {run_filter}")
        # import traceback
        # traceback.print_stack()
    sra_run_ids_filter = set(run_filter) if run_filter is not None else None
    sra_uids = sra_uids_query(query)
    sra_runs = sra_metadata_query(sra_uids, sra_run_ids_filter)
    return sra_runs


def save_sra_metadata_file(dest, sras, non_sras):
    time.sleep(1)
    metadata = list()
    if sras:
        metadata = get_sra_metadata(sras).values()

    sra_meta_keys = ["SRA run accession", "SRA sample accession", "Run type", "SRA read count", "SRA base count", "Average insert size", "Insert size stdev",
                     "Platform type", "Model", "SRA Experiment accession", "SRA Study accession", "Biosample accession", "Bioproject ID",
                     "Bioproject Accession", "Scientific name", "TaxID", "Release date"]
  
    seen_sras = set()
    with open(dest, 'wt') as outf:
        print(f"#{TAB.join(sra_meta_keys)}", file=outf)
        for md in metadata:
            print(TAB.join(str(x) for x in md), file=outf)
            seen_sras.add(md.run_accession)
            
        for ns in non_sras:
            if ns and ns.run_accession not in seen_sras:
                print(TAB.join(str(x) for x in ns), file=outf) 
     

def expand_sra_query(query, run_filter=None):
    # get the list of SRA accessions to be downloaded and metadata if available
    if not query:
        return [], {}
    if type(query) == str:
        # either actual query or file name
        if os.path.exists(query):
            return [], {}
        else:
            # Runs an actual Entrez query and returns the list of SRA accessions
            metadata, _ = sra_query_cached(query, run_filter)

            return list(metadata.keys()), metadata
    elif type(query) == list:
        # list of SRA accessions
        return query, {}
    return query, {}


def find_sra_name(fn):
    return os.path.splitext(os.path.basename(fn))[0].split('_')[0]


def check_for_sra_tools():
    tools = []
    for tool in ['fasterq-dump', 'fastq-dump']:
        if shutil.which(tool):
            tools.append(tool)
    return tools


def download_sra_query(inputs, query_key, sra_dir):
    query = inputs.get(query_key+"_query")
    run_filter = inputs.get(query_key+"_filter")
    is_long_read = query_key == "long_reads"
    sra_runs_list, metadata = expand_sra_query(query, run_filter)
    #print(f"runs: {sra_runs_list}, metadata: {metadata}")
    res = download_sra_runs(sra_runs_list, sra_dir, is_long_read)
    if res != 0:
        return res

    # update sra_reads yaml file with the map from SRA accession
    # to list of fasta files downloaded for this accession.
    existing_fasta_files = glob.glob(os.path.join(sra_dir, '*.fasta'))

    sra_runs_fn = os.path.join(sra_dir, SRA_RUNS_FILE)
    if os.path.exists(sra_runs_fn):
        sra_runs_dict = yaml.safe_load(open(sra_runs_fn, 'r'))
    else:
        sra_runs_dict = {}
    for sra in sra_runs_list:
        sra_runs_dict[sra] = [os.path.abspath(sra_file) for sra_file in existing_fasta_files if find_sra_name(sra_file) == sra]
    if sra_runs_dict:
        Path(sra_dir).mkdir(parents=True, exist_ok=True)
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
        # Convert SraMetadata to list for YAML serialization
        for k, v in metadata.items():
            sra_metadata_dict[k] = list(v) if isinstance(v, SraMetadata) else v
        with open(sra_metadata_fn, 'w') as f:
            yaml.dump(sra_metadata_dict, f)
            f.flush()
        
        sra_queries_fn = os.path.join(sra_dir, SRA_QUERIES_FILE)
        if os.path.exists(sra_queries_fn):
            sra_queries_dict = yaml.safe_load(open(sra_queries_fn, 'r'))
        else:
            sra_queries_dict = {}
        query_hash = hashlib.sha1(query.encode()).hexdigest()
        sra_queries_dict[query_hash] = { "query": query, "runs": sra_runs_list}
        # sra_queries_dict[query] = sra_reads_list
        with open(sra_queries_fn, 'w') as f:
            yaml.dump(sra_queries_dict, f)
            f.flush()
    return 0


def download_sra_runs(sra_runs_list, sra_dir, is_long_read=False):
    read_type = 'long read' if is_long_read else 'short read'
    print(f"RNA-seq {read_type} query returned {len(sra_runs_list)} SRA runs")
    if len(sra_runs_list) == 0:
        return 0
    # check which fasta files are already there and should not be redownloaded
    existing_fasta_files = glob.glob(os.path.join(sra_dir, '*.fasta'))
    # get the list of sra to be downloaded
    existing_sras = set(map(find_sra_name, existing_fasta_files))
    sra_list = list(filter(lambda x : x not in existing_sras, sra_runs_list))
    if len(sra_list) == 0:
        print(f"All SRA runs for {read_type} are already downloaded")
        return 0
    #write the sra_list to be downloaded in file
    Path(sra_dir).mkdir(parents=True, exist_ok=True)
    number_of_files_before = len([entry for entry in os.listdir(sra_dir) if os.path.isfile(os.path.join(sra_dir, entry))])
    tools = check_for_sra_tools()
    for sra in sra_list:
        print(f"Downloading {sra}")
        defline = ">gnl|SRA|$ac.$si.$ri"
        # common arguments for fasterq-dump or fastq-dump
        args = [ sra, '--skip-technical', '--split-files', '-O', sra_dir ]
        args_for_tool = {
            'fasterq-dump' : ['fasterq-dump'] + args + ['--fasta', '--seq-defline', defline, '--threads', '6'],
            'fastq-dump' : ['fastq-dump'] + args + ['--fasta', '0', '--defline-seq', defline] 
        }
        res = 2 # safeguard if tools is empty
        for tool in tools:
            try:
                res = 0
                subprocess.run(args_for_tool[tool], check=True)
            except subprocess.CalledProcessError as err:
                print(f"Error: {tool} failed")
                res = err.returncode
        if res != 0:
            return res
    number_of_files_after = len([entry for entry in os.listdir(sra_dir) if os.path.isfile(os.path.join(sra_dir, entry))])
    print(f"{number_of_files_after - number_of_files_before} files downloaded for {read_type}")
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
    if args.verbosity >= VERBOSITY_VERBOSE:
        print(f"Downloading support data to {args.local_cache}")
    sra_dir = os.path.abspath(os.path.join(args.local_cache, SRA_DOWNLOAD_FOLDER))
    if args.dry_run:
        if args.filename:
            print(f"Downloading SRA")
        return 0
    os.makedirs(args.local_cache, exist_ok=True)
    download_egapx_data(args.local_cache, args.data_version, args.ftp)
    # Set the global cache directory so load_version_map can find the manifest
    global user_cache_dir
    user_cache_dir = args.local_cache
    # Load the version map to populate data_version_cache so that subsequent
    # taxonomy lookups can find the local database instead of calling the API
    load_version_map(args.data_version)
    if args.filename:
        run_inputs = repackage_inputs(yaml.safe_load(open(args.filename, 'r')))
        if args.output:
            run_inputs['output'] = args.output
        print(f"Downloading SRA")
        res, reads_msg = prepare_reads(run_inputs, args.force, for_download=True)
        if not res:
            print('return before downloads')
            print(reads_msg, end='')
            return 1
        inputs = run_inputs['input']
        sra_res = 0
        if check_for_sra_tools():
            for query_key in ['short_reads', 'long_reads']:
                res = download_sra_query(inputs, query_key, sra_dir)
                if res != 0:
                    return res
        elif inputs.get('short_reads_query') or inputs.get('long_reads_query'):
            print("Warning: neither fasterq-dump nor fastq-dump is installed in path")
            sra_res = 1

        # Download BUSCO lineage if needed.
        ensure_busco_lineage_for_inputs(args.local_cache, inputs)
        return sra_res
    return 0


def _cache_target_for_url(url, cache_dir):
    parsed = urlparse(url)
    path = parsed.path.lstrip('/')
    root_prefix = FTP_EGAP_ROOT_PATH.strip('/') + '/'
    if parsed.netloc == FTP_EGAP_SERVER and path.startswith(root_prefix):
        rel = path[len(root_prefix):]
        return os.path.join(cache_dir, rel)
    base_name = os.path.basename(parsed.path.rstrip('/')) or 'downloaded_input'
    digest = hashlib.sha1(url.encode('utf-8')).hexdigest()[:12]
    return os.path.join(cache_dir, 'staging', digest, base_name)


def _find_cached_foreign_path(url, cache_dir):
    """Return cached local path for a foreign URL if present, else empty string."""
    preferred = _cache_target_for_url(url, cache_dir)
    if os.path.exists(preferred):
        return preferred

    # Backward compatibility for files downloaded before cache staging migration.
    parsed = urlparse(url)
    base_name = os.path.basename(parsed.path.rstrip('/')) or 'downloaded_input'
    digest = hashlib.sha1(url.encode('utf-8')).hexdigest()[:12]
    legacy = os.path.join(cache_dir, 'user_inputs', digest, base_name)
    if os.path.exists(legacy):
        return legacy

    return ''


def convert_value_for_cache_download(value, key, strict, cache_dir, download_plan=None, allow_download=True):
    """Convert paths to local-cache targets and collect remote download plan."""
    def _leaf(v, key_name, strict_mode):
        if not isinstance(v, str) or v == '':
            return v

        scheme = _get_uri_scheme(v)
        if allow_download and scheme in {'http', 'https', 'ftp', 's3'}:
            local_target = _cache_target_for_url(v, cache_dir)
            if download_plan is not None:
                download_plan.append({
                    'source': v,
                    'target': local_target,
                    'key': key_name,
                    'scheme': scheme,
                    'foreign': True,
                })
            return local_target

        if scheme:
            # Keep cloud URIs as-is.
            return v

        if not os.path.exists(v):
            if strict_mode:
                raise OSError(ENOENT, f"File for parameter '{key_name}' doesn't exist", v)
            return v
        return os.path.abspath(v)

    return _convert_nested_paths(value, key, strict, _leaf)


def convert_paths_for_cache_download(run_inputs, cache_dir):
    """Convert path-like inputs to local cache paths and return remote download plan."""
    download_plan = []
    _apply_path_converters(
        run_inputs,
        input_converter=lambda v, k, s: convert_value_for_cache_download(
            v, k, s, cache_dir,
            download_plan=download_plan,
            allow_download=True,
        ),
        output_converter=lambda v, k, s: convert_value_for_cache_download(
            v, k, s, cache_dir,
            download_plan=None,
            allow_download=False,
        ),
    )
    return _dedupe_transfer_plan(download_plan)


def ensure_taxonomy_data(cache_dir, use_ftp=False):
    """Download taxonomy subtree required for lineage queries into local cache."""
    taxonomy_db_rel = get_versioned_path("taxonomy", "taxonomy.sqlite3")
    taxonomy_db_local = os.path.join(cache_dir, taxonomy_db_rel)
    if os.path.exists(taxonomy_db_local):
        return

    taxonomy_rel_dir = os.path.dirname(taxonomy_db_rel)
    if verbosity >= VERBOSITY_DEFAULT:
        print(f"Downloading taxonomy data {taxonomy_rel_dir}")
    if use_ftp:
        ftpd = FtpDownloader()
        ftpd.connect(FTP_EGAP_SERVER)
        ftpd.download_ftp_dir(f"{FTP_EGAP_ROOT_PATH}/{taxonomy_rel_dir}", os.path.join(cache_dir, taxonomy_rel_dir))
    else:
        downloader = HttpsDownloader(verbosity=verbosity)
        downloader.download_dir(taxonomy_rel_dir, os.path.join(cache_dir, taxonomy_rel_dir))


def download_needed_data(args):
    """Download only run-required data to local cache for offline execution."""
    if not args.local_cache:
        print("Local cache not set, please use -lc option")
        return 1
    if not args.filename:
        print("Input file must be set for --download-needed mode")
        return 1

    if args.dry_run:
        print("Download-needed dry run")
        return 0

    os.makedirs(args.local_cache, exist_ok=True)

    # Keep manifest/version map behavior compatible with full download mode.
    manifest, effective_data_version = download_data_manifest(args.data_version)
    if args.data_version != effective_data_version:
        print(f"For replicating this run use --data-version {effective_data_version}")
    write_data_manifest_to_cache(args.local_cache, manifest, args.data_version, effective_data_version)
    args.data_version = effective_data_version

    effective_data_version, err_msg = load_version_map(args.data_version)
    if err_msg:
        print(err_msg)
        return 1

    # Taxonomy must be available before expanding parameters that rely on lineage.
    ensure_taxonomy_data(args.local_cache, args.ftp)

    # Activate cache only after taxonomy is present, so downstream lineage calls
    # use local DB and avoid throttled remote fallbacks.
    global user_cache_dir
    user_cache_dir = args.local_cache

    run_inputs = repackage_inputs(yaml.safe_load(open(args.filename, 'r')))
    if args.output:
        run_inputs['output'] = args.output
    if int(args.ortho_taxid) > 0:
        run_inputs['input']['ortho'] = {'taxid': int(args.ortho_taxid)}

    if not expand_and_validate_params(run_inputs):
        return 1

    # Process reads and SRA in the same way as full download mode.
    res, reads_msg = prepare_reads(run_inputs, args.force, for_download=True)
    if not res:
        print(reads_msg, end='')
        return 1

    inputs = run_inputs['input']
    sra_dir = os.path.abspath(os.path.join(args.local_cache, SRA_DOWNLOAD_FOLDER))
    sra_res = 0
    if check_for_sra_tools():
        for query_key in ['short_reads', 'long_reads']:
            res = download_sra_query(inputs, query_key, sra_dir)
            if res != 0:
                return res
    elif inputs.get('short_reads_query') or inputs.get('long_reads_query'):
        print("Warning: neither fasterq-dump nor fastq-dump is installed in path")
        sra_res = 1

    # Handle BUSCO exactly as in full download mode.
    ensure_busco_lineage_for_inputs(args.local_cache, inputs)

    # Convert paths of processed config and cache remote sources into -lc.
    try:
        download_plan = convert_paths_for_cache_download(run_inputs, args.local_cache)
    except OSError as e:
        print(F"{e.strerror}: {e.filename}")
        return 1

    ok, stage_msg = stage_inputs_bulk(download_plan, verbosity=args.verbosity)
    if not ok:
        print(stage_msg)
        return 1
    if args.verbosity >= VERBOSITY_VERBOSE and stage_msg:
        print(stage_msg)

    if args.verbosity > VERBOSITY_QUIET:
        print(reads_msg, end='')
    return sra_res


def main(argv):
    "Main script for EGAPx"
    global verbosity
    # Parse command line
    args = parse_args(argv)

    if args.report_version:
        print(f"EGAPx {get_software_version()}")
        return 0
    verbosity = args.verbosity

    # Make sure than args.data_version and global DATA_VERSION are in sync,
    # preference to explicitely set argument args.data_version
    if not args.data_version:
        args.data_version = DATA_VERSION
    
    if args.download_only or args.download_needed:
        if args.download_needed:
            return download_needed_data(args)
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

    # Preload data version map
    effective_data_version, err_msg = load_version_map(args.data_version)
    if err_msg:
        print(err_msg)
        return 1
    if args.data_version != effective_data_version:
        print(f"For replicating this run use --data-version {effective_data_version}")

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
    if inputs.get('busco_lineage_download'):
        ok, resolved_busco_dir, err = resolve_busco_lineage_download_path(inputs['busco_lineage_download'], busco_lineage)
        if not ok:
            print(f"ERROR: {err}")
            return 1
        inputs['busco_lineage_download'] = resolved_busco_dir
        if args.verbosity >= VERBOSITY_VERBOSE:
            print(f"Using BUSCO lineage from input path: {resolved_busco_dir}")
    elif cache_dir:
        busco_lineage_dir_name = os.path.join(cache_dir, 'busco_downloads', 'lineages', busco_lineage)
        print(f"Checking for BUSCO lineage {busco_lineage} in cache at {busco_lineage_dir_name}")
        if check_busco_lineage(cache_dir, busco_lineage):
            inputs['busco_lineage_download'] = busco_lineage_dir_name
        else:
            print(f"BUSCO lineage {busco_lineage} not found in {cache_dir} or the data is not valid. Please run")
            print(f"  egapx.py -lc {cache_dir} -dl {args.filename}\nto download the lineage")
            return 1

    # Create output directory if needed
    if 'output' in run_inputs and run_inputs['output']: 
        os.makedirs(os.path.join(run_inputs['output'], 'nextflow'), exist_ok=True)

    # Reformat reads into pairs in fromPairs format and add reads_metadata.tsv file
    res, reads_msg = prepare_reads(run_inputs, args.force)
    if not res:
        print(reads_msg, end='')
        return 1

    # Convert inputs to absolute paths and stage foreign inputs.
    try:
        staging_plan = convert_paths_for_executor(
            run_inputs,
            executor=args.executor,
            workdir=workdir,
            staging_dir=args.staging_dir,
        )
    except OSError as e:
        print(F"{e.strerror}: {e.filename}")
        return 1

    if not args.dry_run:
        ok, stage_msg = stage_inputs_bulk(staging_plan, verbosity=args.verbosity)
        if not ok:
            print(stage_msg)
            return 1
        if args.verbosity >= VERBOSITY_VERBOSE and stage_msg:
            print(stage_msg)

    # Add to default task parameters, if input file has some task parameters they will override the default
    task_params = merge_params(task_params, run_inputs)
    ##exit(1)
    task_params['tasks']['cmsearch']['enabled'] = run_inputs['input']['cmsearch']['enabled']
    task_params['tasks']['trnascan']['enabled'] = run_inputs['input'].get('trnascan',dict()).get('enabled', False)

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

    # GP-40510
    # Add commit and branch info to the command-line, so we can filter `nextflow log` to find buildruns of interest.
    in_git_repo = subprocess.run(
            ['git', 'rev-parse', '--is-inside-work-tree'],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        ).returncode == 0
    if in_git_repo:
        nf_cmd += ["--git.commit", subprocess.check_output(['git', 'rev-parse',                 'HEAD'], text=True).strip()]
        nf_cmd += ["--git.branch", subprocess.check_output(['git', 'rev-parse', '--abbrev-ref', 'HEAD'], text=True).strip()]

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
            #print(nf_cmd)
            
            with open(args.output+"/nextflow/nf.config.out", "w") as f:
                config_cmd = ['nextflow', '-C', config_file, 'config' ]
                x = subprocess.run(config_cmd, stdout=f,check=True, text=True, env=env)

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
