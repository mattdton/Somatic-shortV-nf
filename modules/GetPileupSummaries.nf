#!/usr/bin/env nextflow

process GetPileupSummaries_T {

        tag "GetPileupSummaries $bam_id"
        publishDir "${params.outDir}", mode:'copy'

        input:
                path common_biallelic
                path common_biallelic_idx
		tuple val(bam_id) , file(bam_N), file(bam_T)

        output:
                path ("${bam_id}-T.pileups.table")

        shell:
        // Tabulate pileup metrics for inferring contamination - Tumor samples
        '''
        gatk --java-options "-Xmx10g -XX:ParallelGCThreads=2" \
                GetPileupSummaries \
                -I !{bam_T} \
                -V !{common_biallelic} \
                -L !{common_biallelic} \
                -O !{bam_id}-T.pileups.table
        '''
}


process GetPileupSummaries_N {

        tag "GetPileupSummaries $bam_id"
        publishDir "${params.outDir}", mode:'copy'

        input:
                path common_biallelic
                path common_biallelic_idx
		tuple val(bam_id) , file(bam_N), file(bam_T)

        output:
                path ("${bam_id}-N.pileups.table")

        shell:
        // Tabulate pileup metrics for inferring contamination - Normal samples
        '''
        gatk --java-options "-Xmx10g -XX:ParallelGCThreads=2" \
                GetPileupSummaries \
                -I !{bam_N} \
                -V !{common_biallelic} \
                -L !{common_biallelic} \
                -O !{bam_id}-N.pileups.table
        '''

}

