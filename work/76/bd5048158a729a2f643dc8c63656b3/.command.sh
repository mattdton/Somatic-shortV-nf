#!/bin/bash -ue
# GatherVcfs requires intervals in order
ls Patient2-T_Patient2-N.unfiltered.*.vcf.gz   >Patient2_gathered_vcfs_across_subintervals.list

gatk GatherVcfs                 -I  Patient2_gathered_vcfs_across_subintervals.list                 -O  Patient2-T_Patient2-N.unfiltered_unsorted.vcf.gz

#Sort
gatk SortVcf                 -I Patient2-T_Patient2-N.unfiltered_unsorted.vcf.gz                 -O Patient2-T_Patient2-N.unfiltered.vcf.gz
