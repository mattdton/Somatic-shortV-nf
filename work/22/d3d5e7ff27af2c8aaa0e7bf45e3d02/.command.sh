#!/bin/bash -ue
gatk --java-options "-Xmx10g -XX:ParallelGCThreads=2"                 GetPileupSummaries                 -I Patient2-N.coordSorted.dedup.bam                 -V small_exac_common_3.hg38.vcf.gz                 -L small_exac_common_3.hg38.vcf.gz                 -O Patient2-N.pileups.table
