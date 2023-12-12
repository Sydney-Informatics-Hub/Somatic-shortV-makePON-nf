#!/usr/bin/env nextflow

// Enable DSL-2 syntax
nextflow.enable.dsl=2


process Combine_normalCallsUsing_CreateSomaticPanelOfNormals {

	tag "CreateSomaticPanelOfNormals"
        publishDir "${params.outDir}", mode:'copy'
        

        input:
	        path pon_db

        output:
                path 'pon.vcf.gz'
                path 'pon.vcf.gz.tbi'

        script :

        """

        gatk CreateSomaticPanelOfNormals \
                -R ${params.ref}/hs38DH.fasta \
                -V gendb://${pon_db} \
                -O pon.vcf.gz


        """


}
