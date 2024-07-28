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
    # TASK DEFINITIONS
    #
    task BaseRecalibrator{
        input{
            File input_bam
            String recalibration_report_filename
            Array[String] sequence_group_interval
            File dbsnp_vcf
            File dbsnp_vcf_index
            Array[File] known_indels_sites_vcfs
            Array[File] known_indels_sites_indices
            File ref_dict
            File ref_fasta
            File ref_fasta_index
            Int bqsr_scatter
            Int preemptible_tries
            String output_filename

            # Running parameters
            String docker
            String gatk_path
            Int? mem_gb
            Int? disk_space_gb
            Boolean use_ssd = false
        }
        #Float output_bam_size = size(input_cram,"GB")/0.60
        #Int machine_mem_gb = select_first([mem_gb,7])
        #Int command_mem_gb = machine_mem_gb - 1
        Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
        Float dbsnp_size = size(dbsnp_vcf,"GB")
        Int disk_size = ceil(size(input_bam, "GB")/bqsr_scatter) + ref_size +\
                        dbsnp_size + 20

        command{
            ~{gatk_path} --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal\ 
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms5g" \
            BaseRecalibrator\
            -R ~{ref_fasta}\
            -I ~{input_bam}\
            -use-original-qualities \
            --known-sites ~{dbsnp_vcf}\
            --known-sites ~{sep=" -known-sites " known_indels_sites_vcfs}\
            -L -{sep=" -L " sequence_group_interval}\
            -O ~{recalibration_report_filename}
        }
        runtime{
            docker: docker
            preemptible: preemptible_tries
            memory: machine_mem_gb + " GB"
            disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
        }
        output{
            File recalibration_report ="~{recalibration_report_filename}"
        }
    }
    task ApplyBQSR{
        input{
            File input_bam
            String output_bam_basename
            File recalibration_report
            Array[String] sequence_group_interval
            File ref_dict
            File ref_fasta
            File ref_fasta_index
            Int compression_level
            Int bqsr_scatter
            Int preemptible_tries
            String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.19.0"
            # File BaseRecalibrator.output_table
            String output_filename
            String? tmp_dir_name
            File output_table

            # Running parameters
            String gatk_path
            Int? mem_gb
            Int? disk_space_gb
        }
        Float output_bam_size = size(input_cram,"GB")/0.60
        Int machine_mem_gb = select_first([mem_gb,7])
        Int command_mem_gb = machine_mem_gb - 1
        Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
        Int disk_size = ceil(size(input_cram, "GB") + output_bam_size + ref_size) + 20

        command{
            ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G"\
            ApplyBQSR\
            -R ~{ref_fasta}\
            -I ~{input_cram}\
            --bqsr-recal-file BaseRecalibrator.output_table
            --tmp_dir tmp_dir_name
            -O ~{output_filename}.dedup.recal.bam


        }
        runtime{}
        output{}
    }
    
    task CramToBamTask{
        input {
            # Command Parameteres
            File ref_fasta
            File ref_fasta_index
            File ref_dict
            File input_cram
            String sample_name

            # Run parameters
            String docker
            Int? machine_mem_gb
            Int? disk_space_gb
            Boolean use_ssd = false
            Int? preemptible_attempts
            String samtools_path
        }
        Float output_bam_size = size(input_cram,"GB")/0.60
        Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
        Int disk_size = ceil(size(input_cram, "GB") + output_bam_size + ref_size) + 20
        
        command {
            set -e
            set -o pipefail

            ~{samtools_path} view -h -T ~{ref_fasta} ~{input_cram} |
            ~{samtools_path} view -b -o ~{sample_name}.bam -
            ~{samtools_path} index -b ~{sample_name}.bam
            mv ~{sample_name}.bam.bai ~{sample_name}.bai
        }

        runtime {
            docker:docker
            memory: select_first([machine_mem_gb,15])+" GB"
            disks: "local-disk "+select_first([disk_space_gb,disk_size])+if use_ssd then " SSD" else " HDD"
            preemptible: select_first([preemptible_attempts, 3])
        }
        output {
            File output_bam = "~{sample_name}.bam" 
            File output_bai = "~{sample_name}.bai"
        }
    }    

    task HaplotypeCaller {
        input {
            #Command Parameters
            File input_bam
            File input_bam_index
            File interval_list
            String output_filename
            File ref_dict
            File ref_fasta
            File ref_fasta_index
            Float? contamination
            Boolean make_gvcf
            Boolean make_bamout
            Int hc_scatter

            String? gcs_project_for_requester_pays
            
            String gatk_path
            String? java_options

            # Runtime parameters
            String docker
            Int? mem_gb
            Int? disk_space_gb
            Boolean use_ssd = false
            Int? preemptible_attempts
        }    
        
        String java_opt = select_first([java_options, "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"])

        Int machine_mem_gb = select_first([mem_gb,7])
        Int command_mem_gb = machine_mem_gb -1

        Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")
        Int disk_size = ceil(((size(input_bam, "GB") + 30) / hc_scatter) + ref_size) + 20

        String vcf_basename = if make_gvcf then  basename(output_filename, ".gvcf") else basename(output_filename, ".vcf")
        String bamout_arg = if make_bamout then "-bamout ~{vcf_basename}.bamout.bam" else ""

        parameter_meta{
            input_bam: {
                description: "a bam file",
                localization_optional: true
            }    
            input_bam_index: {
                description: "an index file for the bam input",
                localization_optional: true
            }
        }

        command {
            set -e 
            
            ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G ~{java_opt}" \
                HaplotypeCaller \
                -R ~{ref_fasta} \
                -I ~{input_bam} \
                -L ~{interval_list} \
                -O ~{output_filename} \
                -contamination ~{default="0" contamination} \
                -G StandardAnnotation -G StandardHCAnnotation ~{true="-G AS_StandardAnnotation" false="" make_gvcf}\
                -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
                ~{true="-ERC GVCF" false="" make_gvcf}\
                ~{bamout_arg}
            # Cromwell doesn't like optional task outputs, so we have to touch this file.
            touch ~{vcf_basename}.bamout.bam
        }

        runtime {
            docker:docker
            memory: machine_mem_gb + " GB"
            disks: "local-disk "+select_first([disk_space_gb,disk_size]) + if use_ssd then " SSD" else " HDD"
            preemptible: select_first([preemptible_attempts, 3])
        }

        output {
            File output_vcf = "~{output_filename}"
            File output_vcf_index = "~{output_filename}.tbi"
            File bamout = "~{vcf_basename}.bamout.bam"
        }
    }

    # Merge GVCFs
    task MergeGVCFs {
        input{
            #Command parameters
            Array[File] input_vcfs
            Array[File] input_vcfs_indexes
            String output_filename

            String gatk_path

            #Runtime parameters
            String docker
            Int? mem_gb
            Int? disk_space_gb
            Int? preemptible_attempts
        }
            Boolean use_ssd = false
            Int machine_mem_gb = select_first([mem_gb,3])
            Int command_mem_gb = machine_mem_gb - 1
        command{
            set -e 
            ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G" \
            MergeVcfs \
            --INPUT ~{sep=' --INPUT ' input_vcfs}\
            --OUTPUT ~{output_filename}
        }
        runtime{
            docker: docker
            memory: machine_mem_gb + " GB"
            disks: "local-disk " + select_first([disk_space_gb, 100]) + if use_ssd then " SSD" else " HDD"
            preemptible: select_first([preemptible_attempts, 3])
        }
        output{
            File output_vcf = "~{output_filename}"
            File output_vcf_index = "~{output_filename}.tbi"
            
        }
    }

## WORKFLOW DEFINITION
    workflow bqsr_hc_gvcf {
        input {
            File input_bam
            File input_bam_index
            File ref_dict
            File ref_fasta
            File ref_fasta_index
            File known_indels
            File gs_indels
            File dbsnp139
            File scattered_calling_intervals_list

            Boolean make_gvcf = true
            Boolean make_bamout = false
            String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.1.9.0"
            String gatk_path = "/gatk/gatk"
            String gitc_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
            String samtools_path = "samtools"
            String tmp_dir_path = "~{/tmp/}"
        }
        
        Array[File] scattered_calling_intervals = read_lines(scattered_calling_intervals_list)

        Boolean is_cram = sub(basename(input_bam),".*\\.","") == "cram"

        String sample_basename = if is_cram then basename(input_bam,".cram") else basename(input_bam,".bam")
        String vcf_basename = sample_basename
        String output_suffix = if make_gvcf then ".g.vcf.gz" else ".vcf.gz"
        String output_filename = vcf_basename + output_suffix

        Int potential_hc_divisor = length(scattered_calling_intervals) - 20
        Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

        # Check if the input file is CRAM format; If true the convert the CRAM to
        # BAM

        if (is_cram){
            call CramToBamTask{
                input:
                    input_cram = input_bam,
                    sample_name = sample_basename,
                    ref_dict = ref_dict,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    docker = gitc_docker,
                    samtools_path = samtools_path
            }    
        }

        # Call Base Recalibration step
        call BaseRecalibrator{
                input:
                    input_cram = input_bam,
                    output_filename = sample_basename,
                    ref_dict = ref_dict,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    known_indels = known_indels,
                    gs_indels = gs_indels,
                    dbsnp139 = dbsnp139
                }
        call ApplyBQSR{
            input:
                input_cram = input_bam,
                ref_fasta = ref_fasta,
                ref_fasta = ref_fasta_index,
                ref_dict = ref_dict,
                output_table = BaseRecalibrator.output_table,
                output_filename = sample_basename,
                tmp_dir_name = tmp_dir_path,
                gatk_path = gatk_path,

        }
 
        # Call variants in parallel over grouped calling intervals
        #
        scatter (interval_file in scattered_calling_intervals) {
            
            #Generate GVCF by interval 
            call HaplotypeCaller {
                input:
                    input_bam = select_first([CramToBamTask.output_bam,input_bam]),
                    input_bam_index = select_first([CramToBamTask.output_bai, input_bam_index]),
                    interval_list = interval_file,
                    output_filename = output_filename,
                    ref_dict = ref_dict,
                    ref_fasta = ref_fasta,
                    ref_fasta_index = ref_fasta_index,
                    hc_scatter = hc_divisor,
                    make_gvcf = make_gvcf,
                    make_bamout = make_bamout,
                    docker = gatk_docker,
                    gatk_path = gatk_path
            }
        }

        # Merge per-interval GVCFs
        call MergeGVCFs {
            input:
                input_vcfs = HaplotypeCaller.output_vcf,
                input_vcfs_indexes = HaplotypeCaller.output_vcf_index,
                output_filename = output_filename,
                docker = gatk_docker,
                gatk_path = gatk_path
        }

        # Outputs that will be retained when execution is complete
        output {
            File output_vcf = MergeGVCFs.output_vcf
            File output_vcf_index = MergeGVCFs.output_vcf_index
        }
    }


