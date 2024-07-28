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
task annot_vep_snv {
        
        input {
            #Input parameters from command line
            File input_snv_file
            String chrNum
            File gnomad_g_chr_v3
            File gnomad_e_chr_lov
            File gnomad_g_chr_lov
            File dbnsfp_chrNum
            File batch_file
            String? batchNum
             
            #Default parameters
            File ref_fasta
            File ref_dict
            File ref_fasta_index
            File cadd_wgs_snv
            File cadd_wgs_indel
            File exac_pli
            File exac
            File revel_file
            File spliceai_snv
            File spliceai_indel
            File ens_ref_merged
            File ref_fasta
            #String? out_path
            #String? batch_num

            # Running parameters
            String vep_docker

            String file_prefix = basename(input_snv_file,".vcf.gz")
        }

        command <<<

            set -x -e -o pipefail

            mkdir -p output/${batchNum}/${file_prefix}
            cd output/${batchNum}/${file_prefix}

            bcftools view -S ~{batch_file} \
            --force-samples ${input_snv_file} \
            -Ob -o tmpFile.${batchNum}.bcf
            bcftools index tmpFile.${batchNum}.bcf

            bcftools view -r ~{chrNum} tmpFile.${batchNum}.bcf \
            -Ob -o tmpFile.${batchNum}.${chrNum}.bcf
            bcftools index tmpFile.${batchNum}.${chrNum}.bcf
            rm tmpFile.${batchNum}.bcf*

        >>>
        output {
            Array[File] recal_table = glob("output/*")
            #Array[File] g_vcf = glob("output/*.g.vcf.gz")
            #Array[File] g_vcf_index = glob("output/*.g.vcf.gz.tbi")
        }
         runtime {
            docker : vep_docker
            dx_timeout: "48H"
        }
        parameter_meta {
            input_snv_file : {
                description: "Input pVCF chrom files",
                patterns: ["*.vcf.gz"],
                stream: true
            }    
        }
}
    
## WORKFLOW DEFINITION
## To DO: Make HaplotypeCaller using nested scatter and gather 

