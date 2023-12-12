#!/usr/bin/env nextflow

// Enable DSL-2 syntax
nextflow.enable.dsl=2

/*
#Consolidate samples with GenomicsDBImport
#https://github.com/Sydney-Informatics-Hub/Somatic-ShortV/blob/master/gatk4_pon_genomicsdbimport.sh
*/


/*
# This file needs to be created on the fly in bash!! - sample_map_vcf.txt
*/

params.sample_map_vcfs = "Somatic-ShortV/nextflow/make_PON/All_6_make_PON/sample_map_vcf.txt"

process Create_GenomicsDB_from_normalMutect2Calls_GenomicsDBImport {

        tag "GenomicsDBImport"
        publishDir "${params.outDir}", mode:'copy'

        input:        
	        path ('*')

        output:
                path 'pon_db'


        script:

        """


         gatk GenomicsDBImport \
                --sample-name-map $sample_vcfs_input \
                --overwrite-existing-genomicsdb-workspace \
                --genomicsdb-workspace-path pon_db \
                --reader-threads ${task.cpus} \
                --tmp-dir ${temp_dir} \
                --intervals "$base_path/Somatic-ShortV/nextflow/make_PON/All_6_make_PON/100M_primary_interval.list"


        """


}