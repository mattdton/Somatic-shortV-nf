#!/usr/bin/env nextflow



process CalculateContamination {

        tag "CalculateContamination $bam_id"
        publishDir "${params.outDir}", mode:'copy'

        input:
		tuple val(bam_id) , file(bam_N), file(bam_T)
                path pileupsTable_T
                path pileupsTable_N

        output:
                path ("${bam_id}-T_${bam_id}-N_contamination.table") , emit: pair_contamination_table
                path ("${bam_id}-T_segments.table"), emit: tumor_T_segments_table

        shell:
        // Calculate the fraction of reads coming from cross-sample contamination
        
        '''
        
        gatk  --java-options "-Xmx16g" \
                CalculateContamination \
                -I !{bam_id}-T.pileups.table \
                -tumor-segmentation !{bam_id}-T_segments.table \
                -matched !{bam_id}-N.pileups.table \
                -O !{bam_id}-T_!{bam_id}-N_contamination.table
        '''

        }


