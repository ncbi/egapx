

# Eukaryotic Genome Annotation Pipeline - External (EGAPx)

EGAPx is the publicly accessible version of the updated NCBI [Eukaryotic Genome Annotation Pipeline](https://www.ncbi.nlm.nih.gov/genome/annotation_euk/process/). 

EGAPx takes an assembly fasta file, a taxid of the organism, and RNA-seq data. Based on the taxid, EGAPx will pick protein sets, and HMM models. The pipeline runs `miniprot` to align protein sequences, and `STAR` to align RNA-seq, to the assembly. Protein alignments and RNA-seq read alignments are then passed to `Gnomon` for gene prediction. In the first step of `Gnomon`, the short alignments are chained together into putative gene models. In the second step, these predictions are further supplemented by _ab-initio_ predictions based on hmm models. The final annotation for the input assembly is produced as a gff file. 

We currently have protein datasets posted for most vertebrates (mammals, sauropsids, ray-finned fishes), hymenoptera, diptera, lepidoptera and choleoptera, and will be adding datasets for more arthropods, vertebrates, and plants in the next couple of months. Fungi, protists, and nematodes are currently out-of-scope for EGAPx pending additional refinements.

**Warning:**
The current version is an alpha release with limited features and organism scope to collect initial feedback on execution. Outputs are not yet complete and not intended for production use. Please open a GitHub [Issue](https://github.com/ncbi/egapx/issues)  if you encounter any problems with EGAPx. You can also write to cgr@nlm.nih.gov to give us your feedback or if you have any questions.  


**Security Notice:**
EGAPx has dependencies in and outside of its execution path that include several thousand files from the [NCBI C++ toolkit](https://www.ncbi.nlm.nih.gov/toolkit), and more than a million total lines of code. Static Application Security Testing has shown a small number of verified buffer overrun security vulnerabilities. Users should consult with their organizational security team on risk and if there is concern, consider mitigating options like running via VM or cloud instance. 

**License:**
See the EGAPx license [here](https://github.com/ncbi/egapx.git).



## Prerequisites

- Docker or Singularity  
- AWS batch, UGE cluster, or a r6a.4xlarge machine (32 CPUs, 256GB RAM) 
- Nextflow v.23.10.1
- Python v.3.9+ 

Notes:
- General configuration for AWS Batch is described in the Nextflow documentation at https://www.nextflow.io/docs/latest/aws.html
- See Nextflow installation at https://www.nextflow.io/docs/latest/getstarted.html

## The workflow files

- Clone the EGAPx repo:
  ```
  git clone https://github.com/ncbi/egapx.git
  cd egapx
  ```

## Input data format

Input to EGAPx is in the form of a YAML file. 

The following two are the _required_ key-value pairs for the input file:

```
genome: path to assembled genome in FASTA format
taxid: NCBI Taxonomy identifier of the target organism
```

The following are the _optional_ key-value pairs for the input file:

- RNA-seq data. Use one of the following options:
  ```
  reads: [ array of paths to reads FASTA files]
  reads_ids: [ array of SRA run ids ]
  reads_query: query for reads SRA
  ```

- A protein set. A taxid-based protein set will be chosen if no protein set is provided.
  ```
  proteins: path to proteins data in FASTA format. 
  ```

- HMM file used in Gnomon training. A taxid-based HMM will be chosen if no HMM file is provided.
  ```
  hmm: path to HMM file
  ```



## Input example
 
- A test example YAML file `./examples/input_D_farinae_small.yaml` is included in the `egapx` folder. Here, the RNA-seq data is provided as paths to the reads FASTA files. These FASTA files are a sampling of the reads from the complete SRA read files to expedite testing. Currently for paired-end data specified by `reads:`, filenames **must** end in .1 and .2


  ```
  genome: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/809/275/GCF_020809275.1_ASM2080927v1/GCF_020809275.1_ASM2080927v1_genomic.fna.gz
  taxid: 6954
  reads:
    - https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/data/Dermatophagoides_farinae_small/SRR8506572.1
    - https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/data/Dermatophagoides_farinae_small/SRR8506572.2
    - https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/data/Dermatophagoides_farinae_small/SRR9005248.1
    - https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/data/Dermatophagoides_farinae_small/SRR9005248.2
  ```

- To specify an array of NCBI SRA datasets using `reads_ids:`
   ```
   reads_ids:
     - SRR8506572
     - SRR9005248
   ```

- To specify an SRA entrez query using `reads_query:`
    ```
    reads_query: 'txid6954[Organism] AND biomol_transcript[properties] NOT SRS024887[Accession] AND (SRR8506572[Accession] OR SRR9005248[Accession] )'
    ```

  **Note:** Both the above examples `reads_ids` and `reads_query` will have more RNA-seq data than the `input_D_farinae_small.yaml` example. To make sure there aren't a large number of SRA runs, please run your `reads_query` query in NCBI SRA page first. If there are too many SRA runs, then select a few of them and use the `reads_ids` option.   

- First, test EGAPx on the example provided (`input_D_farinae_small.yaml`, a dust mite) to make sure everything works. This example usually runs under 30 minutes depending upon resource availability. There are other examples you can try: `input_C_longicornis.yaml`, a green fly, and `input_Gavia_tellata.yaml`, a bird. These will take close to two hours.  You can prepare your input YAML file following these examples.  

## Run EGAPx

- The `egapx` folder contains the following directories:
    - examples
    - nf
    - test
    - third_party_licenses
    - ui 
     
- The runner script is within the ui directory (`ui/egapx.py`). 

- Create a virtual environment where you can run EGAPx. There is a `requirements.txt` file. PyYAML will be installed in this environment.
  ```
  python -m venv /path/to/new/virtual/environment
  source /path/to/new/virtual/environment/bin/activate
  pip install -r ui/requirements.txt
  ```  


 


- Run EGAPx for the first time to copy the config files so you can edit them:
  ```
  python3 ui/egapx.py ./examples/input_D_farinae_small.yaml  
  ```
  - When you run `egapx.py` for the first time it copies the template config files to the directory `./egapx_config`.
  - You will need to edit these templates to reflect the actual parameters of your setup.
    - For AWS Batch execution you need to edit the value for `process.queue` in `aws.config` file.
    - For local execution on the local machine you don't need to adjust anything.


- Run EGAPx with the following command for real this time.  
  ```
  python3 ui/egapx.py ./examples/input_D_farinae_small.yaml -e aws -w s3://temp_datapath/D_farinae -o example_out
  ```
    
    - use `-e aws` for AWS batch using Docker image
    - use `-e docker` for using Docker image locally
    - use `-e singularity` for using the Singularity image locally
    - type `python3 ui/egapx.py  -h ` for the help menu 

      ```
      $ ./egapx.py  -h

      !!WARNING!!
      This is an alpha release with limited features and organism scope to collect initial feedback on execution. Outputs are not yet complete and not intended for production use.

      usage: egapx.py [-h] [-e EXECUTOR] [-c CONFIG_DIR] [-o OUTPUT] [-w WORKDIR] [-r REPORT] [-n] [-q] [-v] [-fn FUNC_NAME] filename

      Main script for EGAPx

      positional arguments:
        filename              YAML file with input: section with at least genome: and reads: parameters

      optional arguments:
        -h, --help            show this help message and exit
        -e EXECUTOR, --executor EXECUTOR
                              Nextflow executor, one of local, docker, aws. Uses corresponding Nextflow config file
        -c CONFIG_DIR, --config-dir CONFIG_DIR
                              Directory for executor config files, default is ./egapx_config. Can be also set as env EGAPX_CONFIG_DIR
        -o OUTPUT, --output OUTPUT
                              Output path
        -w WORKDIR, --workdir WORKDIR
                              Working directory for cloud executor
        -r REPORT, --report REPORT
                              Report file prefix for report (.report.html) and timeline (.timeline.html) files, default is in output directory
        -n, --dry-run
        -q, --quiet
        -v, --verbose
        -fn FUNC_NAME, --func_name FUNC_NAME
                        func_name
      ```


## Test run

```
$ python3 ui/egapx.py examples/input_D_farinae_small.yaml -e aws -o example_out -w s3://temp_datapath/D_farinae

!!WARNING!!
This is an alpha release with limited features and organism scope to collect initial feedback on execution. Outputs are not yet complete and not intended for production use.

N E X T F L O W  ~  version 23.10.1
Launching `/../home/user/egapx/ui/../nf/ui.nf` [golden_mercator] DSL2 - revision: c134f40af5
in egapx block
executor >  awsbatch (67)
[f5/3007b8] process > egapx:setup_genome:get_genome_info            [100%] 1 of 1 ✔
[32/a1bfa5] process > egapx:setup_proteins:convert_proteins         [100%] 1 of 1 ✔
[96/621c4b] process > egapx:miniprot:run_miniprot                   [100%] 1 of 1 ✔
[6d/766c2f] process > egapx:paf2asn:run_paf2asn                     [100%] 1 of 1 ✔
[56/f1dd6b] process > egapx:best_aligned_prot:run_best_aligned_prot [100%] 1 of 1 ✔
[c1/ccc4a3] process > egapx:align_filter_sa:run_align_filter_sa     [100%] 1 of 1 ✔
[e0/5548d0] process > egapx:run_align_sort                          [100%] 1 of 1 ✔
[a8/456a0e] process > egapx:star_index:build_index                  [100%] 1 of 1 ✔
[d5/6469a6] process > egapx:star_simplified:exec (1)                [100%] 2 of 2 ✔
[64/99ab35] process > egapx:bam_strandedness:exec (2)               [100%] 2 of 2 ✔
[98/a12969] process > egapx:bam_strandedness:merge                  [100%] 1 of 1 ✔
[78/0d7007] process > egapx:bam_bin_and_sort:calc_assembly_sizes    [100%] 1 of 1 ✔
[74/bb014e] process > egapx:bam_bin_and_sort:bam_bin (2)            [100%] 2 of 2 ✔
[39/3cdd00] process > egapx:bam_bin_and_sort:merge_prepare          [100%] 1 of 1 ✔
[01/f64e38] process > egapx:bam_bin_and_sort:merge (1)              [100%] 1 of 1 ✔
[aa/47a002] process > egapx:bam2asn:convert (1)                     [100%] 1 of 1 ✔
[45/6661b3] process > egapx:rnaseq_collapse:generate_jobs           [100%] 1 of 1 ✔
[64/68bc37] process > egapx:rnaseq_collapse:run_rnaseq_collapse (3) [100%] 9 of 9 ✔
[18/bff1ac] process > egapx:rnaseq_collapse:run_gpx_make_outputs    [100%] 1 of 1 ✔
[a4/76a4a5] process > egapx:get_hmm_params:run_get_hmm              [100%] 1 of 1 ✔
[3c/b71c42] process > egapx:chainer:run_align_sort (1)              [100%] 1 of 1 ✔
[e1/340b6d] process > egapx:chainer:generate_jobs                   [100%] 1 of 1 ✔
[c0/477d02] process > egapx:chainer:run_chainer (16)                [100%] 16 of 16 ✔
[9f/27c1c8] process > egapx:chainer:run_gpx_make_outputs            [100%] 1 of 1 ✔
[5c/8f65d0] process > egapx:gnomon_wnode:gpx_qsubmit                [100%] 1 of 1 ✔
[34/6ab0c9] process > egapx:gnomon_wnode:annot (1)                  [100%] 10 of 10 ✔
[a9/e38221] process > egapx:gnomon_wnode:gpx_qdump                  [100%] 1 of 1 ✔
[bc/8ebca4] process > egapx:annot_builder:annot_builder_main        [100%] 1 of 1 ✔
[5f/6b72c0] process > egapx:annot_builder:annot_builder_input       [100%] 1 of 1 ✔
[eb/1ccdd0] process > egapx:annot_builder:annot_builder_run         [100%] 1 of 1 ✔
[4d/6c33db] process > egapx:annotwriter:run_annotwriter             [100%] 1 of 1 ✔
[b6/d73d18] process > export                                        [100%] 1 of 1 ✔
Waiting for file transfers to complete (1 files)
Completed at: 27-Mar-2024 11:43:15
Duration    : 27m 36s
CPU hours   : 4.2
Succeeded   : 67
```
## Output

Look at the output in the out diectory (`example_out`) that was supplied in the command line. The annotation file is called `accept.gff`. 
```
accept.gff
annot_builder_output
nextflow.log
run.report.html
run.timeline.html
run.trace.txt
run_params.yaml
```
The `nextflow.log` is the log file that captures all the process information and their work directories. `run_params.yaml` has all the parameters that was used in the EGAPx run. More information about the process time and resources can be found in the other run* files.  



## Intermediate files

In the above log, each line denotes the process that completed in the workflow. The first column (_e.g._ `[96/621c4b]`) is the subdirectory where the intermediate output files and logs are found for the process in the same line, _i.e._, `egapx:miniprot:run_miniprot`. To see the intermediate files for that process, you can go to the work directory path that you had supplied and traverse to the subdirectory `96/621c4b`: 

```
$ aws s3 ls s3://temp_datapath/D_farinae/96/      
                           PRE 06834b76c8d7ceb8c97d2ccf75cda4/
                           PRE 621c4ba4e6e87a4d869c696fe50034/
$ aws s3 ls s3://temp_datapath/D_farinae/96/621c4ba4e6e87a4d869c696fe50034/
                           PRE output/
2024-03-27 11:19:18          0 
2024-03-27 11:19:28          6 .command.begin
2024-03-27 11:20:24        762 .command.err
2024-03-27 11:20:26        762 .command.log
2024-03-27 11:20:23          0 .command.out
2024-03-27 11:19:18      13103 .command.run
2024-03-27 11:19:18        129 .command.sh
2024-03-27 11:20:24        276 .command.trace
2024-03-27 11:20:25          1 .exitcode
$ aws s3 ls s3://temp_datapath/D_farinae/96/621c4ba4e6e87a4d869c696fe50034/output/
2024-03-27 11:20:24   17127134 aligns.paf
```

## Contact us

Please open a GitHub [Issue](https://github.com/ncbi/egapx/issues)  if you encounter any problems with EGAPx. You can also write to cgr@nlm.nih.gov to give us your feedback or if you have any questions. 
