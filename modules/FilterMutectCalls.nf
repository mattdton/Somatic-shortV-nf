#!/usr/bin/env nextflow



process FilterMutectCalls {

        tag "FilterMutectCalls $bam_id"
        publishDir "${params.outDir}", mode:'copy'



        input :
		tuple val(bam_id) , file(bam_N), file(bam_T)
                path pair_unfiltered_stats
                path pair_contaminationTable
                path segments_T_table
                path pair_unfiltered_vcf
                path pair_unfiltered_vcf_index
                path pair_read_orientation_model
                path(ref)
                

        output :
                path ("${bam_id}-T_${bam_id}-N.filtered.vcf.gz")

        shell:
        // Filter somatic SNVs and indels called by Mutect2
        
        '''
        

        gatk --java-options "-Xmx8g -Xms8g" \
                FilterMutectCalls \
                --reference !{params.ref}/hs38DH.fasta \
                -V !{bam_id}-T_!{bam_id}-N.unfiltered.vcf.gz \
                --stats !{bam_id}-T_!{bam_id}-N.unfiltered.stats \
                --tumor-segmentation !{bam_id}-T_segments.table \
                --contamination-table  !{bam_id}-T_!{bam_id}-N_contamination.table \
                --ob-priors !{bam_id}-T_!{bam_id}-N.read-orientation-model.tar.gz \
                -O !{bam_id}-T_!{bam_id}-N.filtered.vcf.gz



        '''


        }

