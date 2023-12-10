#!/bin/bash -ue
echo "The values: Patient1 Patient1-N.coordSorted.dedup.bam Patient1-T.coordSorted.dedup.bam"

gatk Mutect2              -R /g/data/er01/SIH-HPC-WGS/Reference/hs38DH.fasta              -I Patient1-T.coordSorted.dedup.bam              -I Patient1-N.coordSorted.dedup.bam              -normal Patient1-N              --f1r2-tar-gz Patient1-T_Patient1-N.f1r2.k.tar.gz              -XL chrM              -L /scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-83-Somatic-ShortV/Somatic-shortV-nf/modules/scatter_files/100M_primary_interval_k.list 	     -O Patient1-T_Patient1-N.unfiltered.k.vcf.gz
