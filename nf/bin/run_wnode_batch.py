#!/usr/bin/env python3

import argparse
import subprocess
import shutil
import sys
import os

parser = argparse.ArgumentParser(description="Invoke `command` (must be a worker_node type) on a subset of provided ids.")
parser.add_argument("--ids",                    required=True, help="Seq-ids in scope (all or subset-of) the contents of asn_cache.")
parser.add_argument("--asn-cache",              required=True, help="Passed to gpx_qsubmit and the worker-node.")
parser.add_argument("--batch-num",   type=int,  required=True, help="1..num_batches (inclusive).")
parser.add_argument("--num-batches", type=int,  required=True, help="Total number of batches we are splitting the ids into.")
parser.add_argument("--work-dir",               required=True, help="Temporary output emitted by the workers in the worker-node.")
parser.add_argument("--out-file",               required=True, help="Concatenated output of the worker-node from all workers.")
parser.add_argument("--exclusive", action="store_true",        help="Run in exclusive mode on host.")
parser.add_argument("command", nargs=argparse.REMAINDER,       help="worker_node command followed by command-specific args; common args appended (see script)")
args = parser.parse_args()

assert 1 <= args.num_batches
assert 1 <= args.batch_num <= args.num_batches

os.makedirs(args.work_dir + "/inp/", exist_ok=False)
os.makedirs(args.work_dir + "/out/", exist_ok=False)
jobs_file = args.work_dir + "/inp/jobs.xml"

# generate jobs.xml
subprocess.run(
    ["gpx_qsubmit", "-ids", args.ids, "-nogenbank", "-asn-cache", args.asn_cache, "-o", jobs_file],
    check=True, capture_output=True,
)

# load jobs.xml (one job ber line)
jobs = [line for line in open(jobs_file).read().splitlines() if line and not line.startswith('#')]

# emit sample of jobs corresponding to this batch-id.
with open(args.work_dir + "/inp/jobs_batch.xml", "wt") as f:
    print(
        *(line for i, line in enumerate(jobs) if (i % args.num_batches) + 1 == args.batch_num),
        sep="\n", file=f
    )

# gpx_qdump silently ignores outputs from jobs with job-ids already seen when combining outputs.
# Therefore we can't start assining job-ids from 1 in parallel invocations of wnodes,
# and must ensure that all job-ids are unique (among all invocations of a wnode for a task).
#
# NB: an alternative to this acrobatics is to allow multiple job-ids in gpx_qdump and gpx_make_outputs.
batch_size = -(len(jobs) // -args.num_batches)  # ceildiv
starting_job_id = batch_size * (args.batch_num - 1) + 1

subprocess.run(
    (["flock", "-x", f"/tmp/egapx.{args.command[0]}.lock" ] if args.exclusive else [])
    + args.command
    + [
        "-input-jobs"   , args.work_dir + "/inp/jobs_batch.xml",
        "-start-job-id" , str(starting_job_id),
        "-O"            , args.work_dir + "/out/",
        "-asn-cache"    , args.asn_cache,
        "-nogenbank"
    ], check=True
)

# copy all output to out-file and clean up after ourselves.
with open(args.out_file, "wb") as fout:
    for file in os.listdir(args.work_dir + "/out/"):
        with open(args.work_dir + "/out/" + file, "rb") as f:
            shutil.copyfileobj(f, fout)

shutil.rmtree(args.work_dir)
