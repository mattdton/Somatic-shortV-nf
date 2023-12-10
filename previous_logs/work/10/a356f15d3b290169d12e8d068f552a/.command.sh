#!/bin/bash -ue
gatk IndexFeatureFile                 --input Patient1-T_Patient1-N.filtered.vcf.gz

gatk SelectVariants                 -R /g/data/er01/SIH-HPC-WGS/Reference/hs38DH.fasta                 -V Patient1-T_Patient1-N.filtered.vcf.gz                 --exclude-filtered true                 -O Patient1-T_Patient1-N.filtered_only.vcf.gz
