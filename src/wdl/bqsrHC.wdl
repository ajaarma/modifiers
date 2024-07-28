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
task BaseRecalibrator {
        
        input {
            File input_cram
            String recalibration_report_filename
            File dbsnp_vcf
            File dbsnp_vcf_index
            Array[File] known_indels_sites_vcfs
            Array[File] known_indels_sites_index
            File ref_dict
            File ref_fasta
            File ref_fasta_index
            Int? preemptible_tries

            # Running parameters
            String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.19.0"
            Int? mem_gb
            Int? disk_space_gb
            Boolean use_ssd = false
        }
        Float output_bam_size = size(input_cram,"GB")
        Int machine_mem_gb = select_first([mem_gb,7])
        Int command_mem_gb = machine_mem_gb - 1
        Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
        Float dbsnp_size = size(dbsnp_vcf,"GB")
        Int disk_size = output_bam_size + ref_size + dbsnp_size + 20

        command {
            gatk --java-options "-Xms~{command_mem_gb}G" \
            BaseRecalibrator \
            -R ~{ref_fasta} \
            -I ~{input_cram} \
            --known-sites ~{dbsnp_vcf} \
            --known-sites ~{sep=" -known-sites " known_indels_sites_vcfs} \
            -O ~{recalibration_report_filename}
        }
        runtime {
            docker: gatk_docker
            preemptible: preemptible_tries
            memory: machine_mem_gb + " GB"
            disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
        }
        output {
            File recal_table ="~{recalibration_report_filename}.table"
        }
}
    
task ApplyBQSR {
        input {
            File input_cram
            String output_bam_basename
            File recalibration_report
            File ref_dict
            File ref_fasta
            File ref_fasta_index
            Int? preemptible_tries
            String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.19.0"
            # File BaseRecalibrator.output_table
            Int? mem_gb
            Int additional_disk = 20

        }
        Float output_bam_size = size(input_cram,"GB")
        Int machine_mem_gb = select_first([mem_gb,7])
        Int command_mem_gb = machine_mem_gb - 1
        Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
        Int disk_size = output_bam_size + ref_size + additional_disk

       command {
            gatk --java-options "-Xms~{command_mem_gb}G" \
            ApplyBQSR\
            --create-output-bam-md5 \
            -R ~{ref_fasta} \
            -I ~{input_cram} \
            -O ~{output_bam_basename}.recal.bam \
            -bqsr ~{recalibration_report} \
        }
        runtime {
            docker: gatk_docker
            preemptible: preemptible_tries
            memory: machine_mem_gb+" GB"
            disks: "local-disk " + disk_size + " HDD"
        }
        output {
            File recal_bam = "~{output_bam_basename}.recal.bam"
            File recal_bam_index = "~{output_bam_basename}.recal.bam.bai"
            
        }
}

task HaplotypeCaller {
        input {
            File input_cram
            String output_bam_basename
            File target_intervals
            File ref_dict
            File ref_fasta
            File ref_fasta_index
            File dbsnp_vcf
            String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.19.0"

        }
        #Float output_bam_size = size(input_cram,"GB")
        #Int machine_mem_gb = select_first([mem_gb,7])
        #Int command_mem_gb = machine_mem_gb - 1
        #Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
        #Int disk_size = output_bam_size + ref_size + additional_disk

       command {
            set -x -e -o pipefail

            gatk --java-options "-Xmx2G" \
            HaplotypeCaller\
            -I input_cram \
            -O ~{output_bam_basename}.g.vcf.gz \
            -R ~{ref_fasta} \
            -L ~{target_intervals} \
            -D dbsnp_vcf \
            -G StandardAnnotation -G StandardHCAnnotation \
            -LE true -ERC GVCF -VS LENIENT \
            -native-pair-hmm-threads 2

            rm input_cram
        }
        runtime {
            docker: gatk_docker
            dx_timeout: "48H"
        }
        output {
            File output_gvcf = "~{output_bam_basename}.g.vcf.gz"
            File output_gvcf_index = "~{output_bam_basename}.g.vcf.gz.tbi"
            
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

            String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.9.0"
        }
        
        Boolean is_cram = sub(basename(input_cram),".*\\.","") == "cram"

        String sample_basename = if is_cram then basename(input_cram,".cram") else basename(input_cram,".bam")
        String vcf_basename = sample_basename

        # Call Base Recalibration step

        # To-DO: Loop over each of the samples and then call these tasks
        call BaseRecalibrator {
                input:
                    input_cram = input_cram,
                    recalibration_report_filename = sample_basename,
                    dbsnp_vcf = dbsnp_vcf,
                    dbsnp_vcf_index = dbsnp_vcf_index,
                    known_indels_sites_vcfs = known_indels_sites_vcfs,
                    known_indels_sites_index = known_indels_sites_index,
                    ref_dict = ref_dict,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index
        }
        
        # ApplyBQSR step
        call ApplyBQSR {
            input:
                input_cram = input_cram,
                output_bam_basename = sample_basename,
                recalibration_report = BaseRecalibrator.recal_table,
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index
        }
        
        # HaplotypeCaller steps
        call HaplotypeCaller {
            input:
                input_cram = ApplyBQSR.recal_bam,
                output_bam_basename = sample_basename,
                target_intervals = target_intervals,
                ref_dict = ref_dict,
                ref_fasta = ref_fasta,
                ref_fasta_index = ref_fasta_index,
                dbsnp_vcf = dbsnp_vcf
        }
        
        # Outputs that will be retained when execution is complete
        output {
            #File recal_bam = ApplyBQSR.recal_bam
            File g_vcf = HaplotypeCaller.output_gvcf
            File g_vcf_index = HaplotypeCaller.output_gvcf_index
        }
}


