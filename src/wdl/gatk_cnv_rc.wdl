version 1.0

## This WDL implemets steps to filter rsIDs from imputed genotype data from UKBB
## Requirements: 
##  - input cram (.cram file)
##  - input cram index
##  - reference file
##  - preprocessed interval file
##  - gatk docker
##
## Outputs:
##  - sample specific read count data per target interval (.tsv)
##
# TASK DEFINITIONS
#
task collect_rc {
        
        input {
            #Input parameters from command line
            Array[File] input_cram_list
            Array[File] input_cram_index
            File interval_file
            File ref_file
            File ref_file_index
            File ref_dict

            #Default parameters

            # Running parameters
            String gatk_docker

            #String file_prefix = basename(input_cram,".cram")
        }

        command <<<

            #Commented out because pipe command giving error
            #set -x -e -o pipefail
            #/bin/bash

            #Environment declaration
            mkdir -p output/
            cd output/
          
            #Merge annoation file with samples GT files
            echo "collect read counts"
            for input_cram in ~{sep=" " input_cram_list};do
                sampleID=$( basename $input_cram ".cram")
                echo "$sampleID"

                gatk CollectReadCounts \
                --input ${input_cram} \
                --intervals ~{interval_file} \
                --reference ~{ref_file} \
                --interval-merging-rule OVERLAPPING_ONLY \
                --format TSV \
                --output ${sampleID}.tsv
                
                gzip -v ${sampleID}.tsv
                rm ${input_cram}
            done

        >>>
        output {
            Array[File] out_tsv = glob("output/*.tsv*")
        }
         runtime {
            docker : gatk_docker
            dx_timeout: "48H"
        }
        parameter_meta {
            input_cram_list : {
                description: "Input raw cram file",
                patterns: ["*.cram"]
                #stream: true
            }    
            input_cram_index : {
                description: "Input cram index file",
                patterns: ["*.crai"]
                #stream: true
            }         
        }
}
    
