#!/usr/bin/env nextflow

// Enable DSL-2 syntax
nextflow.enable.dsl=2


process createSampleMapFile_makePON {
    
        tag "create_sample_map_file $bam_id"
        publishDir "${params.outDir}", mode:'copy'


        input:
            path ('*')
        output:
            path ("sample_map_vcf.txt")

'''
# Create an empty text file to store the paths
        output_file="sample_map_vcf.txt"
        > "$output_file"  # Clear contents or create a new file

        # Loop through files matching the pattern and write their paths to the output file
        count=1
        for file in *_gathered_vcfs_across_subintervals_sorted.vcf.gz; do
                # Extract patient name from the file name
                patient_name="Patient$count"

                # Get the full path of the current directory
                current_dir=$(pwd)

                # Write the formatted path to the output file
                echo -e "${patient_name}\t${file}" >> "$output_file"

                ((count++))  # Increment patient count
        done

'''

}