#!/usr/bin/env nextflow

// Enable DSL-2 syntax
nextflow.enable.dsl=2

// Define the process
/// This process runs mutect2 on a tumor/normal sample-pair 

//cpus 1
//container "${params.gatk4__container}"
    

process mutect2 {

        // Unhash container command below and edit name of container
	// if using Docker/Singularity containers
        //container "${params.container}
        
        // where to publish the outputs
        tag "$bam_id $splitIntervalNumber"
        publishDir "$params.outdirA/", mode:'copy'
        
        
        // See: https://www.nextflow.io/docs/latest/process.html#inputs
	/// each input needs to be placed on a new line
        input:
                //path pon_vcf
                //path pon_vcf_index
                tuple val(bam_id) , file(bam_N), file(bam_T)
                each splitIntervalNumber
		
                //path base_path
		//path refdir

        // See: https://www.nextflow.io/docs/latest/process.html#outputs
	// each new output needs to be placed on a new line
        output:
                path ("${bam_id}-T_${bam_id}-N.unfiltered.${splitIntervalNumber}.vcf.gz")
                path ("${bam_id}-T_${bam_id}-N.unfiltered.${splitIntervalNumber}.vcf.gz.stats")
                path ("${bam_id}-T_${bam_id}-N.f1r2.${splitIntervalNumber}.tar.gz")

        // this is the code block 
        script:

        """
        gatk Mutect2 \
             -R $refdir/hs38DH.fasta \
             -I ${bam_N} \
             -I ${bam_T} \
             -normal ${bam_id}-N \
             --f1r2-tar-gz ${bam_id}-T_${bam_id}-N.f1r2.${splitIntervalNumber}.tar.gz \
             -XL chrM \
             -L "$base_path/final_raw_scripts_SomaticShortV_rds/final_scripts_runs_DSL2_DONOTDELETE/100M_primary_interval_${splitIntervalNumber}.list" \
	     -O ${bam_id}-T_${bam_id}-N.unfiltered.${splitIntervalNumber}.vcf.gz

        """




}

