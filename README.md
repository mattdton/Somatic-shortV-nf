# Somatic-shortV-nf

<p align="center">
:wrench: This pipeline is currently under development :wrench:
</p> 

  - [Description](#description)
  - [Diagram](#diagram)
  - [User guide](#user-guide)
      - [Infrastructure usage and
        recommendations](#infrastructure-usage-and-recommendations)
  - [Benchmarking](#benchmarking)
  - [Workflow summaries](#workflow-summaries)
      - [Metadata](#metadata)
      - [Component tools](#component-tools)
  - [Additional notes](#additional-notes)
  - [Help/FAQ/Troubleshooting](#helpfaqtroubleshooting)
  - [Acknowledgements/citations/credits](#acknowledgementscitationscredits)

## Description
Somatic-shortV-nf is a pipeline for identifying somatic short variant (SNPs and indels) events in human Illumina short read whole genome sequence data from tumour and matched normal BAM files. The pipeline follows the [GATK's Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-) workflow. 
The pipeline is written in Nextflow and uses Singularity/Docker to run containerised tools.

There are two main steps to this workflow 
- First step is to generate a large set of candidate somatic variants using the tool  [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360035531132). 
- The next step is to filter the candidate variants to obtain a more confident set of somatic variant calls. 
- The third (optional) step is to annotate the variants with gene level information using an external genomic variant annotations and functional effect prediction tool called [SnpEff](http://pcingola.github.io/SnpEff/).

## Diagram

<p align="center"> 
<img src="./images/Somatic_variant_calling.png" width="60%">
</p> 


## Workflow description 

## User guide 

To run this pipeline, you will need to prepare your input files, reference data, and clone this repository. Before proceeding, ensure Nextflow is installed on the system you're working on. To install Nextflow, see these [instructions](https://www.nextflow.io/docs/latest/getstarted.html#installation). 

### 1. Prepare inputs

To run this pipeline you will need the following inputs: 

* Paired Tumor-Normal paired BAM files
* Corresponding BAM index files  
* Input sample sheet 

This pipeline processes paired BAM files and is capable of processing multiple samples in parallel. BAM files are expected to be coordinate sorted and indexed (see [Fastq-to-BAM](https://github.com/Sydney-Informatics-Hub/Fastq-to-BAM) for an example of a best practice workflow that can generate these files).  

You will need to create a sample sheet with information about the samples you are processing, before running the pipeline. This file must be **comma-separated** and contain a header and one row per sample. Columns should correspond to sampleID, BAM-N file-path, BAM-T file-path: 

```csv
sampleID,bam-N,bam-T 
SAMPLE1,/data/Bams/sample1-N.bam,/data/Bams/sample1-T.bam
SAMPLE2,/data/Bams/sample2-N.bam,/data/Bams/sample2-T.bam|
``````

When you run the pipeline, you will use the mandatory `--input` parameter to specify the location and name of the input file: 

```
--input /path/to/samples.csv
```

### 2. Prepare the reference materials 

To run this pipeline you will need the following reference files:

* Indexed reference genome in FASTA format 

You will need to download and index a copy of the reference genome you would like to use. Reference FASTA files must be accompanied by a .fai index file. If you are working with a species that has a public reference genome, you can download FASTA files from the [Ensembl](https://asia.ensembl.org/info/data/ftp/index.html), [UCSC](https://genome.ucsc.edu/goldenPath/help/ftp.html), or [NCBI](https://www.ncbi.nlm.nih.gov/genome/doc/ftpfaq/) ftp sites. You can use our [IndexReferenceFasta-nf pipeline](https://github.com/Sydney-Informatics-Hub/IndexReferenceFasta-nf) to generate indexes. 

When you run the pipeline, you will use the mandatory `--ref` parameter to specify the location and name of the reference.fasta file: 

```
--ref /path/to/reference.fasta
```


### 3. Clone this repository 

Download the code contained in this repository with: 

```
git clone https://github.com/Sydney-Informatics-Hub/Somatic-shortV-nf
```

This will create a directory with the following structure: 
```
Somatic-shortV-nf/
├── LICENSE
├── README.md
├── config/
├── images/
├── scripts/
├── main.nf
├── modules/
├── PBS_gadi_runPipeline.sh
└── nextflow.config
```
The important features are: 

* **main.nf** contains the main nextflow script that calls all the processes in the workflow.
* **nextflow.config** contains default parameters to use in the pipeline.
* **modules** contains individual process files for each step in the workflow. 
* **config** contains infrastructure-specific config files (this is currently under development)

### 4. Run the pipeline 

The most basic run command for this pipeline is: 

```
nextflow run main.nf --input  sample.csv --ref reference.fasta

```

By default, this will generate `work` directory, `results` output directory and a `runInfo` run metrics directory in the same location you ran the pipeline from. 

To specify additional optional tool-specific parameters, see what flags are supported by running:

```
nextflow run main.nf --help 
```

If for any reason your workflow fails, you are able to resume the workflow from the last successful process with `-resume`. 


### 5. Results 
Once the pipeline is complete, you will find all outputs for each sample in the results directory. **Note** I still need to delete arrange the output files in proper sub-directory structure.

## Infrastructure usage and recommendations

This pipeline has been successfully implemented on [NCI Gadi](https://nci.org.au/our-systems/hpc-systems), [Pawsey Setonix HPC](https://pawsey.org.au/systems/setonix/) and [Pawsey Nimbus cloud](https://pawsey.org.au/systems/nimbus-cloud-service/) using infrastructure-specific configs. These configs can be used to interact with the job scheduler and assign a project code to all task job submissions for billing purposes. You can use the following flags to handle accounting:

* `--whoami` your NCI or Pawsey user name
* `--setonix_account` the Setonix project account you would like to bill service units to
* `--gadi-account` the Gadi project account you would like to bill service units to


### NCI Gadi HPC

Before running the pipeline you will need to load Nextflow and Singularity, both of which are globally installed modules on Gadi. You can do this by running the commands below:

```
module purge
module load nextflow singularity
```

To execute this workflow on NCI Gadi HPC, you will need to specify the following flags to the default run command:

```
nextflow run main.nf --input samples.csv --ref /path/to/ref --gadi-account <account> --whoami <username> -profile gadi
```

Please be aware that as of October 2023, NCI Gadi HPC queues do not have external network access. This means you will not be able to pull the workflow code base or containers if you submit your nextflow run command as a job on any of the standard job queues. NCI currently recommends you run your Nextflow head job either in a GNU screen or tmux session from the login node or submit it as a job to the copyq.

The NCI Gadi config currently runs all tasks on the normal queue. This config uses the `--gadi-account` flag to assign a project code to all task job submissions for billing purposes. The version of Nextflow installed on Gadi has been modified to make it easier to specify resource options for jobs submitted to the cluster. See NCI's [Gadi user guide](https://opus.nci.org.au/display/DAE/Nextflow) for more details.


Please use the PBS script `PBS_gadi_runPipeline.sh` present in the main `Somatic-shortV-nf/` directory for testing the pipeline with test datasets. 


## Benchmarking

Coming soon!


## Workflow summaries
### Metadata

|metadata field     | Somatic-shortV-nf / v1.0     |
|-------------------|:--------------------------------- |
|Version            | 1.0                               |
|Maturity           | under development                 |
|Creators           | Tracy Chew, Cali Willet,Nandan Deshpande                    |
|Source             | NA                                |
|License            | GNU General Public License v3.0   |
|Workflow manager   | NextFlow                          |
|Container          | See component tools               |
|Install method     | NA                                |
|GitHub             | https://github.com/Sydney-Informatics-Hub/Somatic-shortV-nf                            |
|bio.tools          | NA                                |
|BioContainers      | NA                                |
|bioconda           | NA                                | 

### Component tools

To run this pipeline you must have Nextflow and Singularity installed on your machine. All other tools are run using containers. 

|Tool         | Version  |
|-------------|:---------|
|Nextflow     |>=20.07.1 |
|Singularity  |   3.11.3       |
|GATK  |   4.4.0.0       |
|SnpEff       |  5.2   |

## Additional notes
Resources
- [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html)

## Help / FAQ / Troubleshooting
* It is essential that the reference genome you're using contains the same chromosomes, contigs, and scaffolds as the BAM files. To confirm what contigs are included in your indexed BAM file, you can use Samtools idxstats:
```
samtools idxstats input.bam | cut -f 1
```

## License(s)

## Acknowledgements/citations/credits

### Authors 
- Tracy Chew (Sydney Informatics Hub, University of Sydney)
- Cali Willet (Sydney Informatics Hub, University of Sydney)
- Nandan Deshpande (Sydney Informatics Hub, University of Sydney)   

### Acknowledgements 

- This pipeline was built using the [Nextflow DSL2 template](https://github.com/Sydney-Informatics-Hub/Nextflow_DSL2_template).  
- Documentation was created following the [Australian BioCommons documentation guidelines](https://github.com/AustralianBioCommons/doc_guidelines).  

### Cite us to support us! 
Acknowledgements (and co-authorship, where appropriate) are an important way for us to demonstrate the value we bring to your research. Your research outcomes are vital for ongoing funding of the Sydney Informatics Hub and national compute facilities. We suggest including the following acknowledgement in any publications that follow from this work:  

The authors acknowledge the technical assistance provided by the Sydney Informatics Hub, a Core Research Facility of the University of Sydney and the Australian BioCommons which is enabled by NCRIS via Bioplatforms Australia. 

