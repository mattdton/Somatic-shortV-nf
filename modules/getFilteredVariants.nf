#!/usr/bin/env nextflow

process getFilteredVariants {

        tag "getFilteredVariants $bam_id"
        publishDir "${params.outDir}/${bam_id}/getFilteredVariants", mode:'copy'

        input:
                tuple val(bam_id) , file(bam_N), file(bam_T)
                path pair_filtered_vcfs_allvariants
                path(ref)
                
        output:
                path ("${bam_id}-T_${bam_id}-N.filtered_only.vcf.gz")

        shell:
        // Select a subset of variants from a VCF file 

        '''        
        gatk IndexFeatureFile \
                --input !{bam_id}-T_!{bam_id}-N.filtered.vcf.gz

        gatk SelectVariants \
                -R !{params.ref} \
                -V !{bam_id}-T_!{bam_id}-N.filtered.vcf.gz \
                --exclude-filtered true \
                -O !{bam_id}-T_!{bam_id}-N.filtered_only.vcf.gz
        '''
        }
