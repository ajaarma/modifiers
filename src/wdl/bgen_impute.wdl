version 1.0

## This WDL implemets steps to filter rsIDs from imputed genotype data from UKBB
## Requirements: 
##  - raw imputed genotype data (.bgen file)
##  - rsIDlist
##  - chrPosList
##
## Outputs:
##  - chromosome specific bgen file 
##
# TASK DEFINITIONS
#
task bgen_impute {
        
        input {
            #Input parameters from command line
            File input_bgen_file
            File input_bgen_index
            File rs_id_file
            File chr_pos_file
            String chrNum

            #Default parameters

            # Running parameters
            String bgen_docker

            String file_prefix = basename(input_bgen_file,".bgen")
        }

        command <<<

            #Commented out because pipe command giving error
            #set -x -e -o pipefail
            #/bin/bash

            #Environment declaration
            mkdir -p output/
            cd output/
          
            #Merge annoation file with samples GT files
            echo "bgenix output"
            bgenix -g ~{input_bgen_file} -incl-rsids ~{rs_id_file} \
            ~{chr_pos_file} > ~{file_prefix}_~{chrNum}.bgen
            bgenix -g ~{file_prefix}_~{chrNum}.bgen -index

        >>>
        output {
            Array[File] filter_bgen = glob("output/*.bgen*")
        }
         runtime {
            docker : bgen_docker
            dx_timeout: "48H"
        }
        parameter_meta {
            input_bgen_file : {
                description: "Input raw chrom bgren file",
                patterns: ["*.bgen"]
                #stream: true
            }    
            input_bgen_index : {
                description: "Input chrom bgen index file",
                patterns: ["*.bgi"]
                #stream: true
            }         
        }
}
    
