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
task bqsr_hc_all {
        
        input {
            File input_cram
            File input_cram_index
            String recal_report_filename
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

        command {

            gatk --java-options "-Xms4G" \
            BaseRecalibrator \
            -R ~{ref_fasta} \
            -I ~{input_cram} \
            --known-sites ~{dbsnp_vcf} \
            --known-sites ~{sep=" -known-sites " known_indels_sites_vcfs} \
            -O ~{recal_report_filename}.recal.table

            gatk --java-options "-Xms4G" \
            ApplyBQSR \
            -R ~{ref_fasta} \
            -I ~{input_cram} \
            -O ~{recal_report_filename}.recal.bam \
            -bqsr ~{recal_report_filename}.recal.table

            gatk --java-options "-Xmx2G" \
            HaplotypeCaller \
            -I ~{recal_report_filename}.recal.bam \
            -O ~{recal_report_filename}.g.vcf.gz \
            -R ~{ref_fasta} \
            -L ~{target_intervals} \
            -D ~{dbsnp_vcf} \
            -G StandardAnnotation -G StandardHCAnnotation \
            -LE true -ERC GVCF -VS LENIENT \
            -native-pair-hmm-threads 2

            rm ~{recal_report_filename}.recal.bam

        }
        runtime {
            docker: gatk_docker
            dx_instance_type: "mem2_ssd2_v2_x2"
        }
        output {
            File recal_table ="~{recal_report_filename}.recal.table"
            File g_vcf = "~{recal_report_filename}.g.vcf.gz"
            File g_vcf_index = "~{recal_report_filename}.g.vcf.gz.tbi"
        }
}
    

## WORKFLOW DEFINITION
workflow bqsrHC {
        input {
            File input_cram
            File input_cram_index
            File ref_dict
            File ref_fasta
            File ref_fasta_index
            Array[File] known_indels_sites_vcfs
            Array[File] known_indels_sites_index
            File dbsnp_vcf
            File dbsnp_vcf_index

            File target_intervals
            String gatk_docker
        }
        
        Boolean is_cram = sub(basename(input_cram),".*\\.","") == "cram"

        String sample_basename = if is_cram then basename(input_cram,".cram") else basename(input_cram,".bam")
        String vcf_basename = sample_basename

        # Call Base Recalibration step

        # To-DO: Loop over each of the samples and then call these tasks

        call bqsr_hc_all {
                input:
                    input_cram = input_cram,
                    input_cram_index = input_cram_index,
                    recal_report_filename = sample_basename,
                    dbsnp_vcf = dbsnp_vcf,
                    dbsnp_vcf_index = dbsnp_vcf_index,
                    known_indels_sites_vcfs = known_indels_sites_vcfs,
                    known_indels_sites_index = known_indels_sites_index,
                    ref_dict = ref_dict,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    target_intervals  = target_intervals,
                    gatk_docker = gatk_docker
        }
        # Outputs that will be retained when execution is complete
        output {
            #File recal_bam = ApplyBQSR.recal_bam
            File recal_table = bqsr_hc_all.recal_table
            File g_vcf = bqsr_hc_all.g_vcf
            File g_vcf_index = bqsr_hc_all.g_vcf_index
        }
}


