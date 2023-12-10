#!/bin/bash -ue
# For LearnOrientationModel
ls Patient2*f1r2.*.tar.gz > Patient2-T_Patient2-N.unfiltered_f1r2.args


gatk LearnReadOrientationModel --java-options "-Xmx10g"                 --input Patient2-T_Patient2-N.unfiltered_f1r2.args                 -O Patient2-T_Patient2-N.read-orientation-model.tar.gz
