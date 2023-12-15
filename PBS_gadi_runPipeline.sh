#!/bin/bash
  
#PBS -P er01 
#PBS -N Somatic-shortV
#PBS -l walltime=02:00:00
#PBS -l ncpus=1
#PBS -l mem=60GB
#PBS -W umask=022
#PBS -q copyq
#PBS -e Somatic-shortV-nf.e
#PBS -o Somatic-shortV-nf.o
#PBS -l wd
#PBS -l storage=scratch/er01+gdata/er01
#PBS -l jobfs=10GB


#Load singularity and nextflow modules
# See: https://opus.nci.org.au/display/DAE/Nextflow
# See: https://opus.nci.org.au/display/Help/Singularity


module load nextflow/22.04.3
module load singularity



export NXF_SINGULARITY_CACHEDIR=/scratch/$PROJECT/$(whoami)/singularity

# Fill in these variables for your run

samples=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-83-Somatic-ShortV/Somatic-shortV-nf/samples.csv
whoami=npd561
ref=/g/data/er01/SIH-HPC-WGS/Reference
path_to_intervalList=${ref}/BQSR_apply_intervals_13
ponvcf=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-83-Somatic-ShortV/Somatic-shortV-nf/pon.vcf.gz


outDir=results
  


# Run the pipeline (remove annotsv if not needed)
nextflow run main.nf -resume \
        --input ${samples} \
        -profile gadi \
        --whoami ${whoami} --gadi_account $PROJECT \
        --ref ${ref} \
        --ponvcf ${ponvcf} \
        --intervalList_path ${path_to_intervalList} \
        --outDir ${outDir}
        
