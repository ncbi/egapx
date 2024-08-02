# Eukaryotic Genome Annotation Pipeline - External (EGAPx) 

EGAPx is the publicly accessible version of the updated NCBI [Eukaryotic Genome Annotation Pipeline](https://www.ncbi.nlm.nih.gov/genome/annotation_euk/process/). 

EGAPx takes an assembly fasta file, a taxid of the organism, and RNA-seq data. Based on the taxid, EGAPx will pick protein sets and HMM models. The pipeline runs `miniprot` to align protein sequences, and `STAR` to align RNA-seq to the assembly. Protein alignments and RNA-seq read alignments are then passed to `Gnomon` for gene prediction. In the first step of `Gnomon`, the short alignments are chained together into putative gene models. In the second step, these predictions are further supplemented by _ab-initio_ predictions based on HMM models. The final annotation for the input assembly is produced as a `gff` file. 

We currently have protein datasets posted that are suitable for most vertebrates and arthropods:
  - Chordata - Mammalia, Sauropsida, Actinopterygii (ray-finned fishes)
  - Insecta - Hymenoptera, Diptera, Lepidoptera, Coleoptera, Hemiptera 
  - Arthropoda - Arachnida, other Arthropoda

We will be adding datasets for plants and other invertebrates in the next couple of months. Fungi, protists and nematodes are currently out-of-scope for EGAPx pending additional refinements.

We currently have protein datasets posted for most vertebrates (mammals, sauropsids, ray-finned fishes) and arthropods. We will be adding datasets for more arthropods, vertebrates and plants in the next couple of months. Fungi, protists and nematodes are currently out-of-scope for EGAPx pending additional refinements.

**Warning:**
The current version is an alpha release with limited features and organism scope to collect initial feedback on execution. Outputs are not yet complete and not intended for production use. Please open a GitHub [Issue](https://github.com/ncbi/egapx/issues)  if you encounter any problems with EGAPx. You can also write to cgr@nlm.nih.gov to give us your feedback or if you have any questions.  


**Security Notice:**
EGAPx has dependencies in and outside of its execution path that include several thousand files from the [NCBI C++ toolkit](https://www.ncbi.nlm.nih.gov/toolkit), and more than a million total lines of code. Static Application Security Testing has shown a small number of verified buffer overrun security vulnerabilities. Users should consult with their organizational security team on risk and if there is concern, consider mitigating options like running via VM or cloud instance. 

**License:**
See the EGAPx license [here](https://github.com/ncbi/egapx/blob/main/LICENSE).



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

- The following are the _required_ key-value pairs for the input file:

  ```
  genome: path to assembled genome in FASTA format
  taxid: NCBI Taxonomy identifier of the target organism 
  reads: RNA-seq data
  ```
  You can obtain taxid from the [NCBI Taxonomy page](https://www.ncbi.nlm.nih.gov/taxonomy).


  - RNA-seq data can be supplied in any one of the following ways:

    ```
    reads: [ array of paths to reads FASTA or FASTQ files]
    reads: [ array of SRA run IDs ]
    reads: [SRA Study ID]
    reads: SRA query for reads
    ```
  - If you are using your local reads, then the FASTA/FASTQ files should be provided in the following format:
    ```
    reads:
     - path_to_Sample1_R1.gz
     - path_to_Sample1_R2.gz
     - path_to_Sample2_R1.gz
     - path_to_Sample2_R2.gz
    ```

  - If you provide an SRA Study ID, all the SRA run ID's belonging to that Study ID will be included in the EGAPx run.    

- The following are the _optional_ key-value pairs for the input file:  

  - A protein set. A taxid-based protein set will be chosen if no protein set is provided.
    ```
    proteins: path to proteins data in FASTA format. 
    ```

  - HMM file used in Gnomon training. A taxid-based HMM will be chosen if no HMM file is provided.
    ```
    hmm: path to HMM file
    ```



## Input example
 
- A test example YAML file `./examples/input_D_farinae_small.yaml` is included in the `egapx` folder. Here, the RNA-seq data is provided as paths to the reads FASTA files. These FASTA files are a sampling of the reads from the complete SRA read files to expedite testing. 


  ```
  genome: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/809/275/GCF_020809275.1_ASM2080927v1/GCF_020809275.1_ASM2080927v1_genomic.fna.gz
  taxid: 6954
  reads:
    - https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/data/Dermatophagoides_farinae_small/SRR8506572.1
    - https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/data/Dermatophagoides_farinae_small/SRR8506572.2
    - https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/data/Dermatophagoides_farinae_small/SRR9005248.1
    - https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/EGAP/data/Dermatophagoides_farinae_small/SRR9005248.2
  ```

- To specify an array of NCBI SRA datasets:
   ```
   reads:
     - SRR8506572
     - SRR9005248
   ```

- To specify an SRA entrez query:
    ```
    reads: 'txid6954[Organism] AND biomol_transcript[properties] NOT SRS024887[Accession] AND (SRR8506572[Accession] OR SRR9005248[Accession] )'
    ```

  **Note:** Both the above examples will have more RNA-seq data than the `input_D_farinae_small.yaml` example. To make sure the entrez query does not produce a large number of SRA runs, please run it first at the [NCBI SRA page](https://www.ncbi.nlm.nih.gov/sra). If there are too many SRA runs, then select a few of them and list it in the input yaml.   

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
  python3 ui/egapx.py ./examples/input_D_farinae_small.yaml -o example_out 
  ```
  - When you run `egapx.py` for the first time it copies the template config files to the directory `./egapx_config`.
  - You will need to edit these templates to reflect the actual parameters of your setup.
    - For AWS Batch execution, set up AWS Batch Service following advice in the AWS link above. Then edit the value for `process.queue` in `./egapx_config/aws.config` file.
    - For execution on the local machine you don't need to adjust anything.


- Run EGAPx with the following command for real this time.
  - For AWS Batch execution, replace temp_datapath with an existing S3 bucket.
  - For local execution, use a local path for `-w` 
  ```
  python3 ui/egapx.py ./examples/input_D_farinae_small.yaml -e aws -w s3://temp_datapath/D_farinae -o example_out
  ```
    
    - use `-e aws` for AWS batch using Docker image
    - use `-e docker` for using Docker image
    - use `-e singularity` for using the Singularity image
    - use `-e biowulf_cluster` for Biowulf cluster using Singularity image
    - use '-e slurm` for using SLURM in your HPC.
        - Note that for this option, you have to edit `./egapx_config/slurm.config` according to your cluster specifications.
    - type `python3 ui/egapx.py  -h ` for the help menu 

      ```
      $ ui/egapx.py  -h
      
      
      !!WARNING!!
      This is an alpha release with limited features and organism scope to collect initial feedback on execution. Outputs are not yet complete and not intended for production use.

      usage: egapx.py [-h] [-o OUTPUT] [-e EXECUTOR] [-c CONFIG_DIR] [-w WORKDIR] [-r REPORT] [-n] [-st]
                [-so] [-dl] [-lc LOCAL_CACHE] [-q] [-v] [-fn FUNC_NAME]
                [filename]

      Main script for EGAPx

      optional arguments:
        -h, --help            show this help message and exit
        -e EXECUTOR, --executor EXECUTOR
                        Nextflow executor, one of docker, singularity, aws, or local (for NCBI
                        internal use only). Uses corresponding Nextflow config file
        -c CONFIG_DIR, --config-dir CONFIG_DIR
                        Directory for executor config files, default is ./egapx_config. Can be also
                        set as env EGAPX_CONFIG_DIR
        -w WORKDIR, --workdir WORKDIR
                        Working directory for cloud executor
        -r REPORT, --report REPORT
                        Report file prefix for report (.report.html) and timeline (.timeline.html)
                        files, default is in output directory
        -n, --dry-run
        -st, --stub-run
        -so, --summary-only   Print result statistics only if available, do not compute result
        -lc LOCAL_CACHE, --local-cache LOCAL_CACHE
                        Where to store the downloaded files
        -q, --quiet
        -v, --verbose
        -fn FUNC_NAME, --func_name FUNC_NAME
                        func_name

      run:
        filename              YAML file with input: section with at least genome: and reads: parameters
        -o OUTPUT, --output OUTPUT
                        Output path

      download:
        -dl, --download-only  Download external files to local storage, so that future runs can be
                        isolated

      
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
The `nextflow.log` is the log file that captures all the process information and their work directories. `run_params.yaml` has all the parameters that were used in the EGAPx run. More information about the process time and resources can be found in the other run* files.  



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

## Offline mode

If you do not have internet access from your cluster, you can run EGAPx in offline mode. To do this, you would first pull the Singularity image, then download the necessary files from NCBI FTP using `egapx.py` script, and then finally use the path of the downloaded folder in the run command. Here is an example of how to download the files and execute EGAPx in the Biowulf cluster. 


- Download the Singularity image:
```
rm egap*sif
singularity cache clean
singularity pull docker://ncbi/egapx:0.2-alpha
```

- Clone the repo:
```
git clone https://github.com/ncbi/egapx.git
cd egapx
```

- Download EGAPx related files from NCBI:
```
python3 ui/egapx.py -dl -lc ../local_cache
```

- Download SRA reads:
```
prefetch SRR8506572
prefetch SRR9005248
fasterq-dump --skip-technical --threads 6 --split-files --seq-defline ">\$ac.\$si.\$ri" --fasta -O sradir/  ./SRR8506572
fasterq-dump --skip-technical --threads 6 --split-files --seq-defline ">\$ac.\$si.\$ri" --fasta -O sradir/  ./SRR9005248

```
You should see downloaded files inside the 'sradir' folder":
```
ls  sradir/
SRR8506572_1.fasta  SRR8506572_2.fasta  SRR9005248_1.fasta  SRR9005248_2.fasta
```
Now edit the file paths of SRA reads files in `examples/input_D_farinae_small.yaml` to include the above SRA files. 

- Run `egapx.py` first to edit the `biowulf_cluster.config`:
```
ui/egapx.py examples/input_D_farinae_small.yaml -e biowulf_cluster -w dfs_work -o dfs_out -lc ../local_cache
echo "process.container = '/path_to_/egapx_0.2-alpha.sif'" >> egapx_config/biowulf_cluster.config
```

- Run `egapx.py`:
```
ui/egapx.py examples/input_D_farinae_small.yaml -e biowulf_cluster -w dfs_work -o dfs_out -lc ../local_cache

```


## References

Buchfink B, Reuter K, Drost HG. Sensitive protein alignments at tree-of-life scale using DIAMOND. Nat Methods. 2021 Apr;18(4):366-368. doi: 10.1038/s41592-021-01101-x. Epub 2021 Apr 7. PMID: 33828273; PMCID: PMC8026399.

Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. Gigascience. 2021 Feb 16;10(2):giab008. doi: 10.1093/gigascience/giab008. PMID: 33590861; PMCID: PMC7931819.

Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: 10.1093/bioinformatics/bts635. Epub 2012 Oct 25. PMID: 23104886; PMCID: PMC3530905.

Li H. Protein-to-genome alignment with miniprot. Bioinformatics. 2023 Jan 1;39(1):btad014. doi: 10.1093/bioinformatics/btad014. PMID: 36648328; PMCID: PMC9869432.

Shen W, Le S, Li Y, Hu F. SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLoS One. 2016 Oct 5;11(10):e0163962. doi: 10.1371/journal.pone.0163962. PMID: 27706213; PMCID: PMC5051824.



## Contact us

Please open a GitHub [Issue](https://github.com/ncbi/egapx/issues)  if you encounter any problems with EGAPx. You can also write to cgr@nlm.nih.gov to give us your feedback or if you have any questions. 
