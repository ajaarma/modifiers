version 1.0

## This WDL implemets steps to annotate pVCF files from DeepVariant Caller
## Requirements: 
##  - raw unannotated gVCF files
##  - Annotated noGT file
##  - Chromosome Number
##  - batch Number
##  - batch File 
##
## Outputs:
##  - filtered set of 
##  - noGT file
##
# TASK DEFINITIONS
#
task annoMergeGT {
        
        input {
            #Input parameters from command line
            File input_snv_file
            File input_snv_index
            File anno_gt_file
            File anno_gt_index
             
            # Running parameters
            #String bcf_docker

            String file_prefix = basename(input_snv_file,".bcf")
        }

        command <<<

            #Commented out because pipe command giving error

            mkdir -p output/
            cd output/
          
            #Merge annoation file with samples GT files
            echo "Merge annotated GT and Samples"
            bcftools annotate -a ~{anno_gt_file} \
            ~{input_snv_file} -c ID,FILTER,INFO \
            -Ob -o ~{file_prefix}.annot.bcf
            bcftools index ~{file_prefix}.annot.bcf

            rm ~{input_snv_file}
            rm ~{input_snv_index}

        >>>
        output {
            File annot_merge = "output/~{file_prefix}.annot.bcf"
            File annot_merge_index = "output/~{file_prefix}.annot.bcf.csi"

        }
         runtime {
            docker : "dx://project-GV63Xq0Jkqky3PqPj7Jp1VPZ:/resources/soft/bcftools-1.19.tar"
            dx_timeout: "48H"
            dx_instance_type : "mem2_ssd1_v2_x2"
        }
        parameter_meta {
            input_snv_file : {
                description: "Input raw chrom pVCF files",
                patterns: ["*.bcf"]
                #stream: true
            }    
            anno_gt_file : {
                description: "Input annotated GT files",
                patterns: ["*.bcf"]
                #stream: true
            }         
        }
}

