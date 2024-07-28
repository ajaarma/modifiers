version 1.0

## This WDL implemets steps pertaining to GATK pipeline (a) ApplyBQSR (b) HaplotypeCaller
## Requirements: 
## - One analysis-ready BAM files for a singlr samples HC
## - Set of variant calling interval lists for the scatter, provided in a file
##
## Outputs:
## - One GVCF file and index
##
## 
##
##
##
##
##
##
##
##
    # TASK DEFINITIONS
    #
task bqsr_hc_ms {
        
        input {
            Array[File]+ input_cram_list
            #Array[File] input_cram_index
            #String recal_report_filename
            File dbsnp_vcf
            File dbsnp_vcf_index
            Array[File] known_indels_sites_vcfs
            Array[File] known_indels_sites_index
            File ref_dict
            File ref_fasta
            File ref_fasta_index
            File target_intervals

            # Running parameters
            String gatk_docker
        }

        command <<<

            set -x -e -o pipefail 
            mkdir output 
            for input_cram in ~{sep=" " input_cram_list}; do
                
                file_prefix=$( basename $input_cram ".cram")
                
                gatk --java-options "-Xms4G" \
                BaseRecalibrator \
                -R ~{ref_fasta} \
                -I ${input_cram} \
                -L ~{target_intervals} \
                --known-sites ~{dbsnp_vcf} \
                --known-sites ~{sep=" -known-sites " known_indels_sites_vcfs} \
                -O output/${file_prefix}.recal.table

                gatk --java-options "-Xms4G" \
                ApplyBQSR \
                -R ~{ref_fasta} \
                -I ${input_cram} \
                -L ~{target_intervals} \
                -O output/${file_prefix}.recal.bam \
                -bqsr output/${file_prefix}.recal.table

                gatk --java-options "-Xmx4G -XX:ParallelGCThreads=1" \
                HaplotypeCaller \
                -I output/${file_prefix}.recal.bam \
                -O output/${file_prefix}.g.vcf.gz \
                -R ~{ref_fasta} \
                -L ~{target_intervals} \
                -D ~{dbsnp_vcf} \
                -G StandardAnnotation -G StandardHCAnnotation \
                -LE true -ERC GVCF -VS LENIENT \
                -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
                -native-pair-hmm-threads 4

                rm output/${file_prefix}.recal.bam*
                rm ${input_cram}
            done

        >>>
        output {
            Array[File] recal_table = glob("output/*.recal.table")
            Array[File] g_vcf = glob("output/*.g.vcf.gz")
            Array[File] g_vcf_index = glob("output/*.g.vcf.gz.tbi")
        }
         runtime {
            docker : gatk_docker
            dx_instance_type: "mem2_ssd1_v2_x2"
        }
        parameter_meta {
            input_cram_list : {
                description: "Input list of all cram files",
                patterns: ["*.cram"]
            }    
        }
}
    
## WORKFLOW DEFINITION
## To DO: Make HaplotypeCaller using nested scatter and gather 

workflow bqsr_hc_ms_parallel {
        
        input {
            Array[File]+ input_cram_list
            #Array[File] input_cram_index
            #String recal_report_filename
            File dbsnp_vcf
            File dbsnp_vcf_index
            Array[File] known_indels_sites_vcfs
            Array[File] known_indels_sites_index
            File ref_dict
            File ref_fasta
            File ref_fasta_index
            File target_intervals

            # Running parameters
            String gatk_docker
        }

        scatter (for input_cram in input_cram_list){
            call recal_bqsr {
                input : 
                    input_cram = input_cram

                
            }
            scatter (for intervals in target_intervals) {
                
                
                
            }    
            
        }

}
