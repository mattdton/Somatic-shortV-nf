#!/bin/bash -ue
#Combine the stats files across the scattered Mutect2 run

ls Patient1*.stats   >Patient1-T_Patient1-N.unfiltered_stats.args

gatk MergeMutectStats                 --stats Patient1-T_Patient1-N.unfiltered_stats.args                 -O Patient1-T_Patient1-N.unfiltered.stats
