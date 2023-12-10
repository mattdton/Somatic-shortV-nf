#!/bin/bash -ue
gatk IndexFeatureFile                 --input Patient2-T_Patient2-N.filtered.vcf.gz

gatk SelectVariants                 -R /g/data/er01/SIH-HPC-WGS/Reference/hs38DH.fasta                 -V Patient2-T_Patient2-N.filtered.vcf.gz                 --exclude-filtered true                 -O Patient2-T_Patient2-N.filtered_only.vcf.gz
