#!/usr/bin/env nextflow

/// To use DSL-2 will need to include this
nextflow.enable.dsl=2

// =================================================================
// main.nf is the pipeline script for a nextflow pipeline
// ===================================================================


// Import subworkflows to be run in the workflow
include { checkInputs                                                    } from './modules/check_cohort_makePON'
include { mutect2_makePON                                                } from './modules/mutect2_makePON'
include { GatherVcfs_makePON                                             } from './modules/GatherVcfs_makePON'
include { createSampleMapFile_makePON                                    } from './modules/createSampleMapFile_makePON'
include { create_GenomicsDB_makePON                                      } from './modules/create_GenomicsDB_makePON'
include { createSomaticPanelOfNormals_makePON                            } from './modules/createSomaticPanelOfNormals_makePON'



/// Print a header for your pipeline 

log.info """\

      ============================
      ============================
      SOMATIC SHORT V MAKEPoN - NF 
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

 Documentation	@ https://github.com/Sydney-Informatics-Hub/Somatic-shortV-makePON-nf

Cite					@ 10.48546/workflowhub.workflow.431.1 ???

 Log issues @ https://github.com/Sydney-Informatics-Hub/Somatic-shortV-makePON-nf/issues

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
		.splitCsv(header: true, sep:",")
		.map { row -> tuple(row.sampleID, file(row.bam_N), file(row.bam_T))}
	

//params.bams = "/scratch/er01/ndes8648/pipeline_work/nextflow/INFRA-83-Somatic-ShortV/Fastq-to-BAM_BASEDIR/Fastq-to-BAM/Dedup_sort/*-{N,T}.coordSorted.dedup.bam"
//bam_pair_ch=Channel.fromFilePairs( params.bams )

	//Run the processes 
	
  
  mutect2_makePON(bam_pair_ch,params.intervalList)

  GatherVcfs_makePON(mutect2_makePON.out.collect(),bam_pair_ch)

  createSampleMapFile_makePON(GatherVcfs_makePON.out[0].collect())
  
  create_GenomicsDB_makePON(createSampleMapFile_makePON.out,GatherVcfs_makePON.out[0].collect(),GatherVcfs_makePON.out[1].collect())

  createSomaticPanelOfNormals_makePON(create_GenomicsDB_makePON.out)



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
