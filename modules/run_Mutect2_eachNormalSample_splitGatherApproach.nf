#!/usr/bin/env nextflow

// Enable DSL-2 syntax
nextflow.enable.dsl=2


/*
params.bams = "$base_path/bams/M10_reads/*-{Nor,Tum}.recal.bam"
*/

/*
(1)
#https://github.com/Sydney-Informatics-Hub/Somatic-ShortV/blob/master/gatk4_pon.sh

(2)
#https://github.com/Sydney-Informatics-Hub/Somatic-ShortV/blob/master/gatk4_pon_gathervcfs.sh

(3)
#https://github.com/Sydney-Informatics-Hub/Somatic-ShortV/blob/master/gatk4_pon_genomicsdbimport.sh

(4)
#https://github.com/Sydney-Informatics-Hub/Somatic-ShortV/blob/master/gatk4_cohort_pon.sh
*/


/*
#https://github.com/Sydney-Informatics-Hub/Somatic-ShortV/blob/master/README.md

# PON-Step1 - Mutect2: Each Normal Sample
#To create a PoN, call on each normal sample in this mode
#https://github.com/Sydney-Informatics-Hub/Somatic-ShortV/blob/master/gatk4_pon.sh
#To be parallelised

#Scatter and create PoN across genomic intervals
#https://github.com/Sydney-Informatics-Hub/Somatic-ShortV/blob/master/gatk4_pon.sh
#

*/


process run_Mutect2_eachNormalSample_splitGatherApproach  {

        tag "Mutect2 $bam_id $splitIntervalNumber"

        publishDir "${params.outDir}", mode:'copy'

        input:
                tuple val(bam_id) , file(bam_N), file(bam_T)
                each splitIntervalNumber

        output:
                path ("${bam_id}_${splitIntervalNumber}_out.vcf.gz")



        script:

        """

        gatk Mutect2 \
             -R ${params.ref}/hs38DH.fasta \
             -I ${bam_N} \
	     --max-mnp-distance 0 \
             -L ${params.intervalList_path}/100M_primary_interval_${splitIntervalNumber}.list \
             -O ${bam_id}_${splitIntervalNumber}_out.vcf.gz
        """


}









