#!/usr/bin/env nextflow

/// To use DSL-2 will need to include this
nextflow.enable.dsl=2

// =================================================================
// main.nf is the pipeline script for a nextflow pipeline
// ===================================================================


// Import subworkflows to be run in the workflow
include { checkInputs                                } from './modules/check_cohort'
include { mutect2                                    } from './modules/mutect2'
include { GatherVcfs_step                            } from './modules/GatherVcfs_step'
include { MergeMutectStats                           } from './modules/MergeMutectStats'
include { LearnReadOrientationModel                  } from './modules/LearnReadOrientationModel'
include { GetPileupSummaries_T; GetPileupSummaries_N } from './modules/GetPileupSummaries'
include{  CalculateContamination                     } from './modules/CalculateContamination'
include { FilterMutectCalls                          } from './modules/FilterMutectCalls'
include { getFilteredVariants                        } from './modules/getFilteredVariants'
  


/// Print a header for your pipeline 

log.info """\

      ============================
      ============================
          SOMATIC SHORT V - NF 
      ============================
      ============================

 -._    _.--'"`'--._    _.--'"`'--._    _.--'"`'--._    _  
    '-:`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-.  '.` :    
  '.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.:    
  : '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '.  
  '   '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.  
         `-..,..-'       `-..,..-'       `-..,..-'       `       


             ~~~~ Version: 1.0 ~~~~
 

 Created by the Sydney Informatics Hub, University of Sydney

 Documentation	@ https://github.com/Sydney-Informatics-Hub/Somatic-shortV-nf

Cite					@ 10.48546/workflowhub.workflow.431.1 ???

 Log issues @ https://github.com/Sydney-Informatics-Hub/Somatic-shortV-nf/issues

 All the default parameters are set in `nextflow.config`

 =======================================================================================
Workflow run parameters 
=======================================================================================

input       : ${params.input}
outDir      : ${params.outDir}
workDir     : ${workflow.workDir}

=======================================================================================

 """

/// Help function 
// This is an example of how to set out the help function that 
// will be run if run command is incorrect (if set in workflow) 
// or missing/  

def helpMessage() {
    log.info"""
  Usage:  nextflow run main.nf --ref reference.fasta

  Required Arguments:
    --input		  Full path and name of sample input file (tsv format).
	  --ref			  Full path and name of reference genome (fasta format).
	
  Optional Arguments:
    --outDir    Specify name of results directory. 


 HPC accounting arguments:

        --whoami                    HPC user name (Setonix or Gadi HPC)
        --gadi_account              Project accounting code for NCI Gadi (e.g. aa00)
        --setonix_account           Project accounting code for Pawsey Setonix (e.g. name1234)


  """.stripIndent()
}

/// Main workflow structure. 

workflow {

// Show help message if --help is run or if any required params are not 
// provided at runtime

        if ( params.help || params.input == false )
	{   
        // Invoke the help function above and exit
              helpMessage()
              exit 1
	} 

	else 
	{
	
  // Define input bam-pair channel
  // Check inputs file exists
	checkInputs(Channel.fromPath(params.input, checkIfExists: true))
	
  // Original - all files are placed in a specific folder
  //bam_pair_ch=Channel.fromFilePairs( params.bams )

   

  // Split cohort file to collect info for each sample
	bam_pair_ch = checkInputs.out
		.splitCsv(header: true, sep:"\t")
		.map { row -> tuple(row.sampleID, file(row.bam_N), file(row.bam_T))}
	

//params.bams = "/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-83-Somatic-ShortV/Fastq-to-BAM_BASEDIR/Fastq-to-BAM/Dedup_sort/*-{N,T}.coordSorted.dedup.bam"
//bam_pair_ch=Channel.fromFilePairs( params.bams )

	//Run the processes 
	
  // Run mutect2 on a Tumor/Normal sample-pair
  // For initial testing PoN is omitted 
  //mutect2(params.ponvcf,params.ponvcf+'.tbi',bam_pair_ch,intervalList,base_path,refdir_path)
  mutect2(bam_pair_ch,params.intervalList)

  // Gather multiple VCF files from a scatter operation into a single VCF file
	GatherVcfs_step(mutect2.out[0].collect(),bam_pair_ch)

  // Combine the stats files across the scattered Mutect2 run
	MergeMutectStats(mutect2.out[1].collect(),bam_pair_ch)

  // Run the gatk LearnReadOrientationModel 
	LearnReadOrientationModel(mutect2.out[2].collect(),bam_pair_ch)
  

  // Tabulate pileup metrics for inferring contamination - Tumor samples
	//GetPileupSummaries_T(params.ref+'/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz',params.ref+'/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi', bam_pair_ch,LearnReadOrientationModel.out.collect())
  GetPileupSummaries_T(params.ref+'/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz',params.ref+'/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi', bam_pair_ch)


  // Tabulate pileup metrics for inferring contamination - Normal samples
  GetPileupSummaries_N(params.ref+'/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz',params.ref+'/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi', bam_pair_ch)
  
  // Calculate the fraction of reads coming from cross-sample contamination
	CalculateContamination(bam_pair_ch,GetPileupSummaries_T.out.collect(),GetPileupSummaries_N.out.collect())

  // Filter somatic SNVs and indels called by Mutect2
	FilterMutectCalls(bam_pair_ch,MergeMutectStats.out[1].collect(),CalculateContamination.out[0].collect(),CalculateContamination.out[1].collect(),GatherVcfs_step.out[0].collect(),GatherVcfs_step.out[1].collect(),LearnReadOrientationModel.out.collect(),params.ref)	

  // Select a subset of variants from a VCF file 
	getFilteredVariants(bam_pair_ch,FilterMutectCalls.out.collect(),params.ref)

  // Annotate the above subsetted VCF file using snpEff (optional - To be included)
  //annotate_with_snpEff(bam_pair_ch,getFilteredVariants.out)

	}}

workflow.onComplete {
  summary = """
=======================================================================================
Workflow execution summary
=======================================================================================

Duration    : ${workflow.duration}
Success     : ${workflow.success}
workDir     : ${workflow.workDir}
Exit status : ${workflow.exitStatus}
outDir      : ${params.outDir}

=======================================================================================
  """
  println summary

}
