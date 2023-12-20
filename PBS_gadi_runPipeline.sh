#!/bin/bash
  
#PBS -P  
#PBS -N Somatic-shortV
#PBS -l walltime=02:00:00
#PBS -l ncpus=1
#PBS -l mem=20GB
#PBS -W umask=022
#PBS -q copyq
#PBS -e Somatic-shortV-nf.e
#PBS -o Somatic-shortV-nf.o
#PBS -l wd
#PBS -l storage=scratch/[project code]+gdata/[project code]
#PBS -l jobfs=10GB

# Load singularity and nextflow modules
module load nextflow/22.04.3
module load singularity


export NXF_SINGULARITY_CACHEDIR=/scratch/$PROJECT/$(whoami)/singularity

# Fill in these variables for your run

samples=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-83-Somatic-ShortV/Somatic-shortV-nf/samples.csv
ponvcf=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-83-Somatic-ShortV/Somatic-shortV-nf/pon.vcf.gz
ref=/g/data/er01/SIH-HPC-WGS/Reference/hs38DH.fasta
small_exac_common=/g/data/er01/SIH-HPC-WGS/Reference/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz

intervalList_path=/g/data/er01/SIH-HPC-WGS/Reference/BQSR_apply_intervals_13

outDir=results
whoami=npd561  


# Run the pipeline (remove annotsv if not needed)
nextflow run main.nf -resume \
        --input ${samples} \
        -profile gadi \
        --whoami ${whoami} --gadi_account $PROJECT \
        --ref ${ref} --small_exac_common ${small_exac_common}\
        --intervalList_path ${intervalList_path} \
        --ponvcf ${ponvcf} \
        --outDir ${outDir} 