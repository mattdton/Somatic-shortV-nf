#!/bin/bash -ue
# For LearnOrientationModel
ls Patient1*f1r2.*.tar.gz > Patient1-T_Patient1-N.unfiltered_f1r2.args


gatk LearnReadOrientationModel --java-options "-Xmx10g"                 --input Patient1-T_Patient1-N.unfiltered_f1r2.args                 -O Patient1-T_Patient1-N.read-orientation-model.tar.gz
