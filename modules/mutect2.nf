#!/usr/bin/env nextflow

// Enable DSL-2 syntax
nextflow.enable.dsl=2

// Define the process
/// This process runs mutect2 on a tumor/normal sample-pair 



process mutect2 {

        
        // where to publish the outputs
        tag "$bam_id $splitIntervalNumber"
        publishDir "${params.outDir}", mode:'copy'
        
        
        
        input:
                tuple val(bam_id) , file(bam_N), file(bam_T)
                each splitIntervalNumber
		path pon_vcf
                path pon_vcf_index
                
        output:
                path ("${bam_id}-T_${bam_id}-N.unfiltered.${splitIntervalNumber}.vcf.gz") 
                path ("${bam_id}-T_${bam_id}-N.unfiltered.${splitIntervalNumber}.vcf.gz.stats") 
                path ("${bam_id}-T_${bam_id}-N.f1r2.${splitIntervalNumber}.tar.gz")

        
        script:
        // Run mutect2 on a Tumor/Normal sample-pair

        """


        echo "The values: $bam_id $bam_N $bam_T"

        gatk Mutect2 \
             -R ${params.ref}/hs38DH.fasta \
             -I ${bam_T} \
             -I ${bam_N} \
             --panel-of-normals ${pon_vcf} \
             -normal ${bam_id}-N \
             --f1r2-tar-gz ${bam_id}-T_${bam_id}-N.f1r2.${splitIntervalNumber}.tar.gz \
             -XL chrM \
             -L ${params.intervalList_path}/100M_primary_interval_${splitIntervalNumber}.list \
	     -O ${bam_id}-T_${bam_id}-N.unfiltered.${splitIntervalNumber}.vcf.gz

        """




}

