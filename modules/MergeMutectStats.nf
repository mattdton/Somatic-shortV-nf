#!/usr/bin/env nextflow

process MergeMutectStats {

        tag "MergeMutectStats"
        publishDir "${params.outDir}", mode:'copy'

        input:
                path ('*') 
		tuple val(bam_id) , file(bam_N), file(bam_T)
                       
        output:
                path ("${bam_id}-T_${bam_id}-N.unfiltered_stats.args")
                path ("${bam_id}-T_${bam_id}-N.unfiltered.stats")

        shell:
        // Combine the stats files across the scattered Mutect2 run

        """
        ls ${bam_id}*.stats   >${bam_id}-T_${bam_id}-N.unfiltered_stats.args

        gatk MergeMutectStats \
                --stats ${bam_id}-T_${bam_id}-N.unfiltered_stats.args \
                -O ${bam_id}-T_${bam_id}-N.unfiltered.stats

        """
        }


