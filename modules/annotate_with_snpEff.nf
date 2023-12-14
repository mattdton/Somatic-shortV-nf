#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process annotate_with_snpEff {

        tag "annotate_with_snpEff $bam_id"
        publishDir "$params.outdir/Mutect2_filtered/", mode:'copy'

        input:
            tuple val(bam_id) , file(bam_N), file(bam_T)
            path ("${bam_id}-T_${bam_id}-N.filtered_only.vcf.gz")
        
        output:
            path("${bam_id}-T_${bam_id}-N.filtered_only.ann.vcf")


shell:

'''
# https://pcingola.github.io/SnpEff/download/

#snpEff download -dataDir /scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-83-Somatic-ShortV/singularity_cache/ GRCh38.86

snpEff -v \
                -o gatk \
                -stats !{bam_id}-T_!{bam_id}-N.filtered_only.ann.html \
                GRCh38.86 \
                -dataDir /scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-83-Somatic-ShortV/singularity_cache/ \
                !{bam_id}-T_!{bam_id}-N.filtered_only.vcf.gz > !{bam_id}-T_!{bam_id}-N.filtered_only.ann.vcf
'''

}

