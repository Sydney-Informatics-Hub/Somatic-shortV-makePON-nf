#!/usr/bin/env nextflow

// Enable DSL-2 syntax
nextflow.enable.dsl=2


/*
# PON-Step2 - GatherVcfs
# https://github.com/Sydney-Informatics-Hub/Somatic-ShortV/blob/master/gatk4_pon_gathervcfs.sh

TO BE ADDED ---- 
path ('*') from vcf_per_interval_ch.collect()

*/


process GatherVcfs_makePON {

        tag "GatherVcfs_step $bam_id"
        publishDir "${params.outDir}", mode:'copy'


        input:
		path ('*')
                tuple val(bam_id) , file(bam_N), file(bam_T)

        output:
                path ("${bam_id}_gathered_vcfs_across_subintervals_sorted.vcf.gz")
                path ("${bam_id}_gathered_vcfs_across_subintervals_sorted.vcf.gz.tbi")
                path ("${bam_id}_gathered_vcfs_across_subintervals.list")


        script:

        """
        ls ${bam_id}*_out.vcf.gz   >${bam_id}_gathered_vcfs_across_subintervals.list

        # GatherVcfs requires intervals in order, so add chrM using MergeVcfs
        gatk GatherVcfs \
                -I  ${bam_id}_gathered_vcfs_across_subintervals.list \
                -O  ${bam_id}_gathered_vcfs_across_subintervals.vcf.gz

        #Sort
        gatk SortVcf \
                -I ${bam_id}_gathered_vcfs_across_subintervals.vcf.gz \
                -O ${bam_id}_gathered_vcfs_across_subintervals_sorted.vcf.gz
        """

}
