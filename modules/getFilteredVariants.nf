#!/usr/bin/env nextflow
  

nextflow.enable.dsl=2

//container "${params.gatk4__container}"

// Substitute this with a singularity container
//command_path="/scratch/wz54/npd561/installations/snpEff/snpEff"


process getFilteredVariants {

        tag "getFilteredVariants $bam_id"
        publishDir "${params.outDir}", mode:'copy'

        input:
                tuple val(bam_id) , file(bam_N), file(bam_T)
                path pair_filtered_vcfs_allvariants
                path(ref)
                


        output:
                path ("${bam_id}-T_${bam_id}-N.filtered_only.vcf.gz")


        shell:

        '''
        
        gatk IndexFeatureFile \
                --input !{bam_id}-T_!{bam_id}-N.filtered.vcf.gz

        gatk SelectVariants \
                -R !{params.ref}/hs38DH.fasta \
                -V !{bam_id}-T_!{bam_id}-N.filtered.vcf.gz \
                --exclude-filtered true \
                -O !{bam_id}-T_!{bam_id}-N.filtered_only.vcf.gz

        '''

        }
