#!/bin/bash
  
#PBS -P er01 
#PBS -N Somatic-shortV
#PBS -l walltime=01:00:00
#PBS -l ncpus=
#PBS -l mem=80GB
#PBS -W umask=022
#PBS -q copyq
#PBS -e germlineStructuralV-nf.e
#PBS -o germlineStructuralV-nf.o
#PBS -l wd
#PBS -l storage=scratch/er01+gdata/er01

#Load singularity and nextflow modules
# See: https://opus.nci.org.au/display/DAE/Nextflow
# See: https://opus.nci.org.au/display/Help/Singularity
module load java
#module load nextflow/21.04.1
module load nextflow/22.04.3
module load singularity
#module load gatk/4.1.8.1
module load gatk/4.1.4.0


whoami=er01

# Run the pipeline (remove annotsv if not needed)
nextflow run main.nf 
        --input ${samples} -profile gadi \
        --refdir ${ref} \
        --whoami ${whoami} --gadi_account $PROJECT