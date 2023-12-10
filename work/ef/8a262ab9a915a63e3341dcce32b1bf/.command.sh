#!/bin/bash -ue
gatk  --java-options "-Xmx16g"                 CalculateContamination                 -I Patient1-T.pileups.table                 -tumor-segmentation Patient1-T_segments.table                 -matched Patient1-N.pileups.table                 -O Patient1-T_Patient1-N_contamination.table
