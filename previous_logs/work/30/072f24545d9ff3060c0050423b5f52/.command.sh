#!/bin/bash -ue
# GatherVcfs requires intervals in order
ls Patient1-T_Patient1-N.unfiltered.*.vcf.gz   >Patient1_gathered_vcfs_across_subintervals.list

gatk GatherVcfs                 -I  Patient1_gathered_vcfs_across_subintervals.list                 -O  Patient1-T_Patient1-N.unfiltered_unsorted.vcf.gz

#Sort
gatk SortVcf                 -I Patient1-T_Patient1-N.unfiltered_unsorted.vcf.gz                 -O Patient1-T_Patient1-N.unfiltered.vcf.gz
