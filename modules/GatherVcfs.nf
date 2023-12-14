#!/usr/bin/env nextflow



process GatherVcfs {

        tag "GatherVcfs $bam_id"
        publishDir "${params.outDir}", mode:'copy'


        input:
                path ('*')
                tuple val(bam_id) , file(bam_N), file(bam_T)
                

        output:
                path ("${bam_id}-T_${bam_id}-N.unfiltered.vcf.gz") 
                path ("${bam_id}-T_${bam_id}-N.unfiltered.vcf.gz.tbi")

        
        script:
        // Gather multiple VCF files from a scatter operation into a single VCF file

        """
        

        # GatherVcfs requires intervals in order
        ls ${bam_id}-T_${bam_id}-N.unfiltered.*.vcf.gz   >${bam_id}_gathered_vcfs_across_subintervals.list
 
        gatk GatherVcfs \
                -I  ${bam_id}_gathered_vcfs_across_subintervals.list \
                -O  ${bam_id}-T_${bam_id}-N.unfiltered_unsorted.vcf.gz

        #Sort
        gatk SortVcf \
                -I ${bam_id}-T_${bam_id}-N.unfiltered_unsorted.vcf.gz \
                -O ${bam_id}-T_${bam_id}-N.unfiltered.vcf.gz

        """

}

