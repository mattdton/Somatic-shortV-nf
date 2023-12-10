#!/bin/bash -ue
gatk  --java-options "-Xmx16g"                 CalculateContamination                 -I Patient2-T.pileups.table                 -tumor-segmentation Patient2-T_segments.table                 -matched Patient2-N.pileups.table                 -O Patient2-T_Patient2-N_contamination.table
