#!/bin/bash

#PBS -P er01
#PBS -N pullContainers
#PBS -l walltime=02:00:00
#PBS -l ncpus=1
#PBS -l mem=190GB
#PBS -W umask=022
#PBS -q copyq
#PBS -l wd
#PBS -l storage=scratch/er01+gdata/er01

module load singularity

# TODO create lint test to check this script and nextflow.config have same containers

# specify singularity cache dir (consistent with gadi.config)
export SINGULARITY_CACHE_DIR=/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-83-Somatic-ShortV/singularity_cache

# pull containers
singularity pull --dir $SINGULARITY_CACHE_DIR docker://quay.io/biocontainers/gatk4:4.1.4.1--0
