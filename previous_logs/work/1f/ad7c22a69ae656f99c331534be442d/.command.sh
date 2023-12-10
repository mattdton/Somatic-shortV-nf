#!/bin/bash -ue
#Combine the stats files across the scattered Mutect2 run

ls Patient2*.stats   >Patient2-T_Patient2-N.unfiltered_stats.args

gatk MergeMutectStats                 --stats Patient2-T_Patient2-N.unfiltered_stats.args                 -O Patient2-T_Patient2-N.unfiltered.stats
