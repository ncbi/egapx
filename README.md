# Eukaryotic Genome Annotation Pipeline - External (EGAPx) 

EGAPx is the publicly accessible version of the updated NCBI [Eukaryotic Genome Annotation Pipeline](https://www.ncbi.nlm.nih.gov/refseq/annotation_euk/process/). 

EGAPx takes an assembly FASTA file, a taxid of the organism, and RNA-seq data. Based on the taxid, EGAPx will pick protein sets and HMM models. The pipeline runs `miniprot` to align protein sequences, and `STAR` to align RNA-seq to the assembly. Protein alignments and RNA-seq read alignments are then passed to `Gnomon` for gene prediction. In the first step of `Gnomon`, the short alignments are chained together into putative gene models. In the second step, these predictions are further supplemented by _ab-initio_ predictions based on HMM models. Functional annotation is added to the final structural annotation set based on the type and quality of the model and orthology information. The final annotation for the input assembly is produced as a `gff` file. 

We currently have protein datasets posted that are suitable for most vertebrates, arthropods, and some plants:
  - Chordata - Mammalia, Sauropsida, Actinopterygii (ray-finned fishes), other Vertebrates
  - Insecta - Hymenoptera, Diptera, Lepidoptera, Coleoptera, Hemiptera 
  - Arthropoda - Arachnida, other Arthropoda

  - Monocots - Lilipopsida
  - Eudicots - Asterids, Rosids, Fabids, Caryophyllales
  

Fungi, protists and nematodes are currently out-of-scope for EGAPx pending additional refinements.

**Submitting to GenBank:**
If you’d like to be an early tester as we refine the output and workflow for submitting EGAPx annotation to GenBank, please contact us at cgr@nlm.nih.gov.

**Warning:**
The current version is an early release and still under active development to add features and refine outputs. The workflow for GenBank submission is still under development. Please open a GitHub [Issue](https://github.com/ncbi/egapx/issues) if you encounter any problems with EGAPx. You can also write to cgr@nlm.nih.gov to give us your feedback or if you have any questions.  


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
  - The assembled genome should be screened for contamination prior to running EGAPx. See the NCBI [Foreign Contamination Screen](https://github.com/ncbi/fcs) for a fast, user-friendly contamination screening tool. 

  - You can obtain taxid from the [NCBI Taxonomy page](https://www.ncbi.nlm.nih.gov/taxonomy).


  - RNA-seq data can be supplied in any one of the following ways:

    ```
    reads: [ array of paths to reads FASTA or FASTQ files]
    reads: [ array of SRA run IDs ]
    reads: [SRA Study ID]
    reads: SRA query for reads
    ```

- The following are the _optional_ key-value pairs for the input file. The default taxid-based settings (i.e. omitting these parameters) are recommended for most use cases:  

  - A protein set. A taxid-based protein set will be chosen if no protein set is provided. This should only be needed for annotation of obscure organisms or those with little RNAseq data available.
    ```
    proteins: path to proteins data in FASTA format. 
    ```

  - HMM file used in Gnomon training. A taxid-based HMM will be chosen if no HMM file is provided.
    ```
    hmm: path to HMM file
    ```

- The following are _optional_ metadata configuration parameters (not critical for testing EGAPx alpha, will later be used for GenBank submissions):
  - Annotation provider. The main contact for this genome assembly.
    ```
    annotation_provider: GenBank submitter 
    ```
  - Annotation name prefix. GenBank assembly accession version. Uniquely identifies annotation products originating from the same annotation run. The resulting annotation name is `<annotation_name_prefix>-GB_YYYY_MM_DD`. If the GCA acc.ver if not known, do not include this parameter. The annotation name will default to `GB_YYYY_MM_DD`.
    ```
    annotation_name_prefix: GCA_#########.1 
    ```
  - Locus tag prefix. One to 9-letter prefix to use for naming genes on this genome assembly. If an official locus tag prefix was already reserved from an INSDC organization (GenBank, ENA or DDBJ) for the given BioSample and BioProject pair, provide here. Otherwise, provide a string of your choice. If no value is provided, the prefix 'egapxtmp' will be used.
    ```
    locus_tag_prefix: egapxtmp 
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


- If you are using your local reads, then the FASTA/FASTQ files can be provided using the format below. For proper specification of paired-end read files, the filenames must have a shared prefix prior to an underscore character, and the prefix is not shared by any other library:
   ```
   reads:
     - path/to/se1_reads.fq      # path to single-end reads
     - path/to/se2_reads.fq
     - path/to/pe1_reads_R1.fq   # path to paired-end R1 reads
     - path/to/pe1_reads_R2.fq   # path to paired-end R2 reads
     - path/to/pe2_reads_R1.fq
     - path/to/pe2_reads_R2.fq
    ```
    
    Alternatively, you can explicitly set the names and paths to reads sets by following the format below. Here the filenames for the reads can be anything, but the set names for each set has to be unique. 
    ```
    reads:
     - - single_end_library_name1   # set name
       - - path/to/se1_reads.fq     # file name for single-end reads
     - - single_end_library_name2
       - - path/to/se2_reads.fq
     - - paired_end_library_name1   # set name  
       - - path/to/pe1_reads_R1.fq  # file name for paired-end R1 reads
         - path/to/pe1_reads_R2.fq  # file name for paied-end R2 reads
     - - paired_end_library_name2
       - - path/to/pe2_reads_R1.fq
         - path/to/pe2_reads_R2.fq
    ```
    
    There is one other option on how to provide RNA-seq data. If you have a large number of local RNA-seq data, you can list them in a file with a set name and a filepath in each line (see `examples/input_D_farinae_small_reads.txt`). Then you can read that file from the input yaml (see `examples/input_D_farinae_small_readlist.yaml`). 
     
- To specify an array of NCBI SRA datasets:
   ```
   reads:
     - SRR8506572
     - SRR9005248
   ```
   - If you provide an SRA Study ID, all the SRA run ID's belonging to that Study ID will be included in the EGAPx run.    

- To specify an SRA entrez query:
    ```
    reads: 'txid6954[Organism] AND biomol_transcript[properties] NOT SRS024887[Accession] AND (SRR8506572[Accession] OR SRR9005248[Accession] )'
    ```

  **Note:** Both the above examples using SRA reads or SRA entrez query will have more RNA-seq data than the `input_D_farinae_small.yaml` example. To make sure the entrez query does not produce a large number of SRA runs, please run it first at the [NCBI SRA page](https://www.ncbi.nlm.nih.gov/sra). If there are too many SRA runs, then select a few of them and list it in the input yaml.   

- First, test EGAPx on the example provided (`input_D_farinae_small.yaml`, a dust mite) to make sure everything works. This example usually runs under 30 minutes depending upon resource availability. There are other examples you can try: `input_C_longicornis.yaml`, a green fly, and `input_Gavia_stellata.yaml`, a bird. These will take close to two hours.  You can prepare your input YAML file following these examples.  

## Run EGAPx

- The `egapx` folder contains the following directories:
    - examples
    - nf
    - third_party_licenses
    - ui 
     
- The runner script is within the ui directory (`ui/egapx.py`). 

- Create a virtual environment where you can run EGAPx. There is a `requirements.txt` file. PyYAML will be installed in this environment.
  ```
  python -m venv /path/to/new/virtual/environment
  source /path/to/new/virtual/environment/bin/activate
  pip install -r requirements.txt
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
    - use `-e slurm` for using SLURM in your HPC.
        - Note that for this option, you have to edit `./egapx_config/slurm.config` according to your cluster specifications.
    - type `python3 ui/egapx.py  -h ` for the help menu 

      ```
      $ ui/egapx.py  -help
      usage: egapx.py [-h] [-o OUTPUT] [-e EXECUTOR] [-c CONFIG_DIR] [-w WORKDIR]
                      [-r REPORT] [-n] [-st] [-so] [-ot ORTHO_TAXID] [-dl]
                      [-lc LOCAL_CACHE] [-q] [-v] [-V] [-fn FUNC_NAME]
                      [filename]

      Main script for EGAPx

      optional arguments:
        -h, --help            show this help message and exit
        -e EXECUTOR, --executor EXECUTOR
                              Nextflow executor, one of docker, singularity, aws,
                              or local (for NCBI internal use only). Uses
                              corresponding Nextflow config file
        -c CONFIG_DIR, --config-dir CONFIG_DIR
                              Directory for executor config files, default is
                              ./egapx_config. Can be also set as env
                              EGAPX_CONFIG_DIR
        -w WORKDIR, --workdir WORKDIR
                              Working directory for cloud executor
        -r REPORT, --report REPORT
                              Report file prefix for report (.report.html) and
                              timeline (.timeline.html) files, default is in output
                              directory
        -n, --dry-run
        -st, --stub-run
        -so, --summary-only   Print result statistics only if available, do not
                              compute result
        -ot ORTHO_TAXID, --ortho-taxid ORTHO_TAXID
                              Taxid of reference data for orthology tasks
        -lc LOCAL_CACHE, --local-cache LOCAL_CACHE
                              Where to store the downloaded files
        -q, --quiet
        -v, --verbose
        -V, --version         Report software version
        -fn FUNC_NAME, --func_name FUNC_NAME
                              func_name

      run:
        filename              YAML file with input: section with at least genome:
                              and reads: parameters
        -o OUTPUT, --output OUTPUT
                        Output path

      download:
        -dl, --download-only  Download external files to local storage, so that
                              future runs can be isolated
      ```


## Test run

```
$ python3 ui/egapx.py examples/input_D_farinae_small.yaml -e aws -o example_out -w s3://temp_datapath/D_farinae

!!WARNING!!
This is an alpha release with limited features and organism scope to collect initial feedback on execution. Outputs are not yet complete and not intended for production use.

N E X T F L O W  ~  version 23.10.1
Launching `/../home/user/egapx/ui/../nf/ui.nf` [golden_mercator] DSL2 - revision: c134f40af5
in egapx block
executor >  awsbatch (83)
[41/69fb92] process > egapx:setup_genome:get_genome_info                                                                  [100%] 1 of 1 ✔
[12/af924a] process > egapx:setup_proteins:convert_proteins                                                               [100%] 1 of 1 ✔
[26/661e33] process > egapx:target_proteins_plane:miniprot:split_proteins                                                 [100%] 1 of 1 ✔
[86/68836c] process > egapx:target_proteins_plane:miniprot:run_miniprot (1)                                               [100%] 1 of 1 ✔
[f1/2d07a3] process > egapx:target_proteins_plane:paf2asn:run_paf2asn (1)                                                 [100%] 1 of 1 ✔
[05/33457c] process > egapx:target_proteins_plane:best_aligned_prot:run_best_aligned_prot                                 [100%] 1 of 1 ✔
[41/455b4f] process > egapx:target_proteins_plane:align_filter_sa:run_align_filter_sa                                     [100%] 1 of 1 ✔
[c9/4627b4] process > egapx:target_proteins_plane:align_sort_sa:run_align_sort                                            [100%] 1 of 1 ✔
[9b/0b248b] process > egapx:rnaseq_short_plane:star_index:build_index                                                     [100%] 1 of 1 ✔
[79/799e31] process > egapx:rnaseq_short_plane:star:run_star (1)                                                          [100%] 2 of 2 ✔
[01/af1f68] process > egapx:rnaseq_short_plane:bam_strandedness:rnaseq_divide_by_strandedness                             [100%] 1 of 1 ✔
[65/4107dc] process > egapx:rnaseq_short_plane:bam_bin_and_sort:calc_assembly_sizes                                       [100%] 1 of 1 ✔
[5d/c69fbf] process > egapx:rnaseq_short_plane:bam_bin_and_sort:bam_bin (2)                                               [100%] 2 of 2 ✔
[c1/707e59] process > egapx:rnaseq_short_plane:bam_bin_and_sort:merge_prepare                                             [100%] 1 of 1 ✔
[e3/bba172] process > egapx:rnaseq_short_plane:bam_bin_and_sort:merge (1)                                                 [100%] 1 of 1 ✔
[2b/7c7b6a] process > egapx:rnaseq_short_plane:bam2asn:convert (1)                                                        [100%] 1 of 1 ✔
[23/3a9fba] process > egapx:rnaseq_short_plane:rnaseq_collapse:generate_jobs                                              [100%] 1 of 1 ✔
[b8/994db8] process > egapx:rnaseq_short_plane:rnaseq_collapse:run_rnaseq_collapse (8)                                    [100%] 9 of 9 ✔
[da/f769f6] process > egapx:rnaseq_short_plane:rnaseq_collapse:run_gpx_make_outputs                                       [100%] 1 of 1 ✔
[af/c32ba6] process > egapx:gnomon_plane:chainer:run_align_sort (1)                                                       [100%] 1 of 1 ✔
[7f/bed27d] process > egapx:gnomon_plane:chainer:generate_jobs                                                            [100%] 1 of 1 ✔
[4a/cdb342] process > egapx:gnomon_plane:chainer:run_chainer (7)                                                          [100%] 16 of 16 ✔
[7c/b687bb] process > egapx:gnomon_plane:chainer:run_gpx_make_outputs                                                     [100%] 1 of 1 ✔
[62/e78572] process > egapx:gnomon_plane:gnomon_wnode:gpx_qsubmit                                                         [100%] 1 of 1 ✔
[62/8445b3] process > egapx:gnomon_plane:gnomon_wnode:annot (1)                                                           [100%] 10 of 10 ✔
[57/589794] process > egapx:gnomon_plane:gnomon_wnode:gpx_qdump                                                           [100%] 1 of 1 ✔
[7b/020592] process > egapx:annot_proc_plane:fetch_swiss_prot_asn                                                         [100%] 1 of 1 ✔
[70/34b131] process > egapx:annot_proc_plane:get_swiss_prot_ids                                                           [100%] 1 of 1 ✔
[7d/16a826] process > egapx:annot_proc_plane:prot_gnomon_prepare:prot_gnomon_prepare_p                                    [100%] 1 of 1 ✔
[a3/a6a568] process > egapx:annot_proc_plane:diamond_worker:run_diamond_egap                                              [100%] 1 of 1 ✔
[97/e54b4a] process > egapx:annot_proc_plane:best_protein_hits:run_protein_filter_replacement                             [100%] 1 of 1 ✔
[e3/32a317] process > egapx:annot_proc_plane:gnomon_biotype:run_gnomon_biotype                                            [100%] 1 of 1 ✔
[89/56953c] process > egapx:annot_proc_plane:annot_builder:annot_builder_main                                             [100%] 1 of 1 ✔
[7c/28df80] process > egapx:annot_proc_plane:annot_builder:annot_builder_input                                            [100%] 1 of 1 ✔
[19/781bc2] process > egapx:annot_proc_plane:annot_builder:annot_builder_run                                              [100%] 1 of 1 ✔
[f5/1140c6] process > egapx:annot_proc_plane:print_fake_lxr_data                                                          [100%] 1 of 1 ✔
[94/0ee74c] process > egapx:annot_proc_plane:orthology_plane:fetch_ortholog_references                                    [100%] 1 of 1 ✔
[f3/053877] process > egapx:annot_proc_plane:orthology_plane:setup_ext_genome:get_genome_info                             [100%] 1 of 1 ✔
[bd/5ededd] process > egapx:annot_proc_plane:orthology_plane:setup_ext_proteins:convert_proteins                          [100%] 1 of 1 ✔
[7d/fa5f13] process > egapx:annot_proc_plane:orthology_plane:get_prot_ref_ids                                             [100%] 1 of 1 ✔
[82/8018fb] process > egapx:annot_proc_plane:orthology_plane:extract_products_from_models:run_extract_products_from_mo... [100%] 1 of 1 ✔
[ce/22bdea] process > egapx:annot_proc_plane:orthology_plane:diamond_orthology:run_diamond_egap                           [100%] 1 of 1 ✔
[ed/0d0cdd] process > egapx:annot_proc_plane:orthology_plane:find_orthologs:run_find_orthologs                            [100%] 1 of 1 ✔
[56/48bd29] process > egapx:annot_proc_plane:locus_track:run_locus_track                                                  [100%] 1 of 1 ✔
[95/4ad706] process > egapx:annot_proc_plane:locus_link:run_locus_link                                                    [100%] 1 of 1 ✔
[1e/a66cb3] process > egapx:annot_proc_plane:final_asn_markup:final_asn                                                   [100%] 1 of 1 ✔
[f2/391794] process > egapx:annot_proc_plane:annotwriter:run_annotwriter                                                  [100%] 1 of 1 ✔
[4e/6fccc1] process > egapx:convert_annotations:run_converter                                                             [100%] 1 of 1 ✔
[8d/e3225f] process > export                                                                                              [100%] 1 of 1 ✔
Completed at: 30-Oct-2024 11:46:09
Duration    : 53m 9s
CPU hours   : 7.0
Succeeded   : 83


Statistics for example_out/complete.genomic.gff
CDS          33203
exon         35007
gene         8828
lnc_RNA      566
mRNA         8407
pseudogene   6
transcript   4
```

## Output

Look at the output in the out diectory (`example_out`) that was supplied in the command line. The annotation file is called `complete.genomic.gff`. 
```
annot_builder_output
annotated_genome.asn
annotation_data.cmt
complete.cds.fna
complete.genomic.fna
complete.genomic.gff
complete.genomic.gtf
complete.proteins.faa
complete.transcripts.fna
nextflow.log
resume.sh
run.report.html
run.timeline.html
run.trace.txt
run_params.yaml
stats
validated
```
Description of the outputs:
* `complete.genomic.gff`: final annotation set in GFF3 format.
* `complete.genomic.gtf`: final annotation set in GTF format.
* `complete.genomic.fna`: full genome sequences set in FASTA format.
* `complete.genomic.gtf`: final annotation set in gtf format.
* `complete.cds.fna`: annotated Coding DNA Sequences (CDS) in FASTA format.
* `complete.transcripts.fna`: annotated transcripts in FASTA format (includes UTRs).
* `complete.proteins.faa`: annotated protein products in FASTA format.
* `annotated_genome.asn`: final annotation set in ASN1 format.

Description of the logs and miscellaneous outputs:
* `annot_builder_output/accept.ftable_annot`: intermediate file with accepted annotation models called by Gnomon.
* `annotation_data.cmt`: annotation structured comment file. Used for submission to GenBank.
* `nextflow.log`: main Nextflow log that captures all the process information and their work directories.
* `resume.sh`: Nextflow command for resuming a run from the last successful task.
* `run.report.html`: Nextflow rendered HTML execution report containing run summary, resource usage, and tasks execution.
* `run.timeline.html`: Nextflow rendered HTML timeline for all processes executed in the EGAPx pipeline.
* `run.trace.txt`: Nextflow execution tracing file that contains information about each EGAPx process including runtime and CPU usage.
* `run_params.yaml`: YAML file containing parameters used for the EGAPx run
* `stats`: directory containing features statistics for the final annotation set
* `validated`: directory containing validation warnings and errors for annotated features. Used for submission to GenBank.

## Interpreting Output

`stats/feature_counts.xml` contains summary counts of features by model prediction categories determined by Gnomon.

**NOTE** not all categories are expected to have counts data (e.g. model RefSeq, fully supported, ab initio)

Genes with `major correction` are likely protein-coding genes with frameshifts and/or internal stops. These models include "LOW QUALITY PROTEIN" in the protein FASTA title, are marked up with exception=low-quality sequence region on the mRNA and CDS features, and the annotation is adjusted to meet GenBank criteria (frameshifts are compensated for by 1-2 bp microintrons in the mRNA and CDS features, and internal stops have a transl_except to translate the codon as X instead of a stop). For RefSeq, we set a threshold of no more than 10% of protein-coding genes with major corrections to release the annotation. We recommend users polish assembly sequences if the rate is higher than 10%.

Counts of protein-coding genes should be considered versus similar species. Low counts may result from insufficient supporting evidence (e.g. low RNAseq coverage or an unusual organism compared to the available protein data). High counts may indicate genome fragmentation, uncollapsed haplotypic duplication, or noise from genes annotated on transposons.

`stats/feature_stats.xml` contains summary statistics on transcript counts per gene, exon counts per transcript, and the counts and length distributions of features by sub-type.

## Intermediate files

In the nextflow log, you can find work directory paths for each job. You can go to that path, and look for the output files and command logs. For example, to see the files generated during run_miniprot job, run the following command to get the directory path, and list the files within that directory.

```
grep run_miniprot example_out/nextflow.log| grep COMPLETED

aws s3 ls s3://temp_datapath/D_farinae/86/68836c310a571e6752a33a221d1962/
                           PRE output/
2024-10-30 10:54:36          0 
2024-10-30 10:59:04          6 .command.begin
2024-10-30 10:59:33        780 .command.err
2024-10-30 10:59:35        780 .command.log
2024-10-30 10:59:32          0 .command.out
2024-10-30 10:54:36      13013 .command.run
2024-10-30 10:54:36        139 .command.sh
2024-10-30 10:59:33        277 .command.trace
2024-10-30 10:59:34          1 .exitcode

aws s3 ls s3://ncbi-egapx-expires/work/D_farinae/86/68836c310a571e6752a33a221d1962/output/
2024-10-30 10:59:34   26539116 1.paf
```

## Offline mode

If you do not have internet access from your cluster, you can run EGAPx in offline mode. To do this, you would first pull the Singularity image, then download the necessary files from NCBI FTP using `egapx.py` script, and then finally use the path of the downloaded folder in the run command. Here is an example of how to download the files and execute EGAPx in the Biowulf cluster. 


- Download the Singularity image:
```
rm egap*sif
singularity cache clean
singularity pull docker://ncbi/egapx:0.3.2-alpha
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
echo "process.container = '/path_to_/egapx_0.3.2-alpha.sif'" >> egapx_config/biowulf_cluster.config
```

- Run `egapx.py`:
```
ui/egapx.py examples/input_D_farinae_small.yaml -e biowulf_cluster -w dfs_work -o dfs_out -lc ../local_cache

```

## Modifying default parameters

The default task parameter values are listed in the file `ui/assets/default_task_params.yaml`. If there are cases where you need to change some task parameters from the default values, you can add those to the input yaml file.  For example, if you're using RNA-seq from species besides the one being annotated, you can relax the alignment criteria by setting the following parameters in your input yaml:

```
tasks:
  rnaseq_collapse:
    rnaseq_collapse: -high-identity 0.8
  convert_from_bam:
    sam2asn: -filter 'pct_identity_gap >= 85'
  star_wnode:
    star_wnode: -pct-identity 85
```



## References

Buchfink B, Reuter K, Drost HG. Sensitive protein alignments at tree-of-life scale using DIAMOND. Nat Methods. 2021 Apr;18(4):366-368. doi: 10.1038/s41592-021-01101-x. Epub 2021 Apr 7. PMID: 33828273; PMCID: PMC8026399.

Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. Gigascience. 2021 Feb 16;10(2):giab008. doi: 10.1093/gigascience/giab008. PMID: 33590861; PMCID: PMC7931819.

Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: 10.1093/bioinformatics/bts635. Epub 2012 Oct 25. PMID: 23104886; PMCID: PMC3530905.

Li H. Protein-to-genome alignment with miniprot. Bioinformatics. 2023 Jan 1;39(1):btad014. doi: 10.1093/bioinformatics/btad014. PMID: 36648328; PMCID: PMC9869432.

Shen W, Le S, Li Y, Hu F. SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLoS One. 2016 Oct 5;11(10):e0163962. doi: 10.1371/journal.pone.0163962. PMID: 27706213; PMCID: PMC5051824.



## Contact us

Please open a GitHub [Issue](https://github.com/ncbi/egapx/issues)  if you encounter any problems with EGAPx. You can also write to cgr@nlm.nih.gov to give us your feedback or if you have any questions. 
