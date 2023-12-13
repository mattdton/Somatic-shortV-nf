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


#module load nextflow/21.04.1
module load java
module load nextflow/22.04.3
module load singularity
module load gatk/4.1.8.1
#module load gatk/4.1.4.0


export NXF_SINGULARITY_CACHEDIR=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-83-Somatic-ShortV/singularity_cache

# Fill in these variables for your run
ref=/g/data/er01/SIH-HPC-WGS/Reference
samples=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-83-Somatic-ShortV/Somatic-shortV-nf/sample.tsv
#annotsvDir= #full path to directory housing Annotations_Human directory
#annotsvMode= #both|full|split. see annotation mode in https://github.com/lgmgeo/AnnotSV/blob/master/README.AnnotSV_3.3.4.pdf
whoami=npd561
path_to_intervalList=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-83-Somatic-ShortV/Somatic-shortV-nf/modules/scatter_files
outDir=results
ponvcf=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-83-Somatic-ShortV/Somatic-shortV-nf/pon.vcf.gz
  


# Run the pipeline (remove annotsv if not needed)
nextflow run main.nf -resume \
        --input ${samples} \
        -profile gadi \
        --whoami ${whoami} --gadi_account $PROJECT \
        --ref ${ref} \
        --ponvcf ${ponvcf} \
        --intervalList_path ${path_to_intervalList} \
        --outDir ${outDir}
        
