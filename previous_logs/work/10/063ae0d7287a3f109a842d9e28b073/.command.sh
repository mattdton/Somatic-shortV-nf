#!/bin/bash -ue
echo "The values: Patient2 Patient2-N.coordSorted.dedup.bam Patient2-T.coordSorted.dedup.bam"

gatk Mutect2              -R /g/data/er01/SIH-HPC-WGS/Reference/hs38DH.fasta              -I Patient2-T.coordSorted.dedup.bam              -I Patient2-N.coordSorted.dedup.bam              -normal Patient2-N              --f1r2-tar-gz Patient2-T_Patient2-N.f1r2.l.tar.gz              -XL chrM              -L /scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-83-Somatic-ShortV/Somatic-shortV-nf/modules/scatter_files/100M_primary_interval_l.list 	     -O Patient2-T_Patient2-N.unfiltered.l.vcf.gz
