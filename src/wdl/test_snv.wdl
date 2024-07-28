version 1.0

## This WDL implemets steps to annotate pVCF files from DeepVariant Caller
## Requirements: 
##  - One analysis-ready gVCF files
##  - Chromosome Number
##  - batch Number
##  - batch File 
##
## Outputs:
##  - batch wise bcf file
##  - noGT file
##
# TASK DEFINITIONS
#
task test_vep_snv {
        
        input {
            #Input parameters from command line
            File input_snv_file
            String chrNum
            String batchNum
            File batch_file
             
            #Default parameters
            File ref_fasta
            File ref_dict
            File ref_fasta_index
            File exon_region

            # Running parameters
            String vep_docker

            String file_prefix = basename(input_snv_file,".vcf.gz")
        }

        command <<<

            #Commented out because pipe command giving error
            #set -x -e -o pipefail
            #/bin/bash

            #Environment declaration
            #export PATH=/opt/vep/bcftools-1.9/bin:$PATH
            #export BCFTOOLS_PLUGINS=/opt/vep/bcftools-1.9/plugins
            #export PATH=/opt/vep/custom:$PATH

            #mkdir -p outputs/~{batchNum}/~{file_prefix}
            #cd outputs/~{batchNum}/~{file_prefix}
            mkdir -p output/
            cd output/
           
            #Extract sample set
            echo "Extract ~{batchNum} samples"
            bcftools view -S ~{batch_file} \
            --force-samples ~{input_snv_file} \
            -Oz -o ~{file_prefix}.~{batchNum}.bcf
            bcftools index ~{file_prefix}.~{batchNum}.bcf

            #Normalize bcfs
            echo "Normalization of variants; Split multi allelic sites"
            bcftools norm -m - -f ~{ref_fasta} ~{file_prefix}.~{batchNum}.bcf | \
            bcftools view -e 'AC==0' -Ob -o ~{file_prefix}.~{batchNum}.norm.bcf
            bcftools index ~{file_prefix}.~{batchNum}.norm.bcf

            if ![ -s ~{file_prefix}.~{batchNum}.norm.bcf ]; then
                echo -e "'ERROR: ~{file_prefix}.~{batch_file}.norm.bcf ' doesn't exist"
                exit 1
            else 
                echo "Normalization step completed"
                echo "Removing File: ~{file_prefix}.~{batchNum}.bcf*"
                rm ~{file_prefix}.~{batchNum}.bcf*
            fi
            
            #Extract exonic regions
            echo "Extract variants for only exonic regions"
            bcftools view -R ~{exon_region} ~{file_prefix}.~{batchNum}.norm.bcf \
            -Ob -o ~{file_prefix}.~{batchNum}.norm.exon.bcf

            bcftools view ~{file_prefix}.~{batchNum}.norm.exon.bcf |\
            bcftools sort -m 1G -T ~/output/ | bcftools norm -D |\
            bcftools view -Ob -o ~{file_prefix}.~{batchNum}.norm.exon.sorted.bcf

            mv ~{file_prefix}.~{batchNum}.norm.exon.sorted.bcf \
                ~{file_prefix}.~{batchNum}.norm.exon.bcf
            bcftools index ~{file_prefix}.~{batchNum}.norm.exon.bcf

            if ![ -s ~{file_prefix}.~{batchNum}.norm.exon.bcf ]; then
                echo -e "'ERROR: ~{file_prefix}.~{batchNum}.norm.exon.bcf ' doesn't exist"
                exit 1
            else 
                echo "Normalization step completed"
                echo "Removing File: ~{file_prefix}.~{batchNum}.norm.bcf*"
                rm ~{file_prefix}.~{batchNum}.norm.bcf*
            fi

            #Extract only genotype
            bcftools view -G ~{file_prefix}.~{batchNum}.norm.exon.bcf \
            -Oz -o ~{file_prefix}.~{batchNum}.norm.exon.noGT.vcf.gz
            bcftools index ~{file_prefix}.~{batchNum}.norm.exon.noGT.vcf.gz
            
            #Annotate gnomAD -genome,exome v2.1,1, genome-v3.0
            
        >>>
        output {
            Array[File] norm_file = glob("output/*.exon.bcf*")
            Array[File] g_vcf = glob("output/*.vcf.gz*")
            #Array[File] g_vcf_index = glob("output/*.g.vcf.gz.tbi")
        }
         runtime {
            docker : vep_docker
            dx_timeout: "48H"
        }
        parameter_meta {
            input_snv_file : {
                description: "Input pVCF chrom files",
                patterns: ["*.vcf.gz"]
                #stream: true
            }    
        }
}
    
## WORKFLOW DEFINITION
## To DO: Make HaplotypeCaller using nested scatter and gather 

