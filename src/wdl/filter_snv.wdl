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
task filter_snv {
        
        input {
            #Input parameters from command line
            File input_snv_file
            File input_snv_index
            File anno_gt_file
            File anno_gt_index
            String chrNum
            String batchNum
            File exon_region
            Array[File] batch_file_list

             
            #Default parameters

            # Running parameters
            String vep_docker

            String file_prefix = basename(input_snv_file,".bcf")
        }

        command <<<

            #Commented out because pipe command giving error
            #set -x -e -o pipefail
            #/bin/bash

            #Environment declaration
            export PATH=/opt/vep/bcftools-1.9/bin:$PATH
            export BCFTOOLS_PLUGINS=/opt/vep/bcftools-1.9/plugins
            export PATH=/opt/vep/custom:$PATH

            #mkdir -p outputs/~{batchNum}/~{file_prefix}
            #cd outputs/~{batchNum}/~{file_prefix}
            mkdir -p output/
            cd output/
          
            #Merge annoation file with samples GT files
            echo "Merge annotated GT and Samples ~{batchNum} samples"
            bcftools annotate -a ~{anno_gt_file} \
            ~{input_snv_file} -c ID,FILTER,INFO \
            -Ob -o ~{file_prefix}.annot.bcf
            bcftools index ~{file_prefix}.annot.bcf

            rm ~{input_snv_file}
            rm ~{input_snv_index}

            #Extract exonic region
            echo "Extract exonic regions"
            bcftools view -R ~{exon_region} \
            ~{file_prefix}.annot.bcf -Ob -o \
            ~{file_prefix}.annot.exon.bcf

            #Sorted exon coordinates
            bcftools view ~{file_prefix}.annot.exon.bcf | \
            bcftools sort -m 6G -T /home/dnanexus/output | \
            bcftools norm -D | bcftools view -Ob -o \
            ~{file_prefix}.annot.exon.sorted.bcf

            mv ~{file_prefix}.annot.exon.sorted.bcf \
            ~{file_prefix}.annot.exon.bcf
            bcftools index ~{file_prefix}.annot.exon.bcf
                
            #Filter against gnomAD MAF <=0.01 or 1%
            if [ $chrNum != "chrY"]; then

                #Filter against gnomAD MAF <=0.01 or 1%
                echo "Filter against gnomAD MAF <=0.01"
                bcftools view -Ob -e 'GNOMADg4_AF>0.01' \
                ~{file_prefix}.annot.exon.bcf | \
                bcftools view -Ob -e 'GNOMADe4_AF>0.01' > \
                ~{file_prefix}.annot.exon.0_01.bcf
            else
                echo "Filter against gnomAD MAF <=0.01 for chrY"
                bcftools view -Ob -e 'GNOMADe4_AF>0.01' \
                ~{file_prefix}.annot.exon.bcf > \
                ~{file_prefix}.annot.exon.0_01.bcf
            fi
 
            bcftools index ~{file_prefix}.annot.exon.0_01.bcf

            rm ~{file_prefix}.annot.exon.bcf*

            #Filter by Impact
            echo "Filter by Impact"
            bcftools view ~{file_prefix}.annot.exon.0_01.bcf | \
            grep -E '^#|^MT|MODERATE|HIGH|splice_donor_5th_base_variant|splice_donor_0th_base_variant|splice_donor_region_variant|splice_polypyrimidine_tract_variant|HGMD|CLNSIG' | \
            bcftools view -Ob -o ~{file_prefix}.annot.exon.0_01.imp.bcf
            bcftools index ~{file_prefix}.annot.exon.0_01.imp.bcf

            rm ~{file_prefix}.annot.exon.0_01.bcf*

            #Convert VCF 2 TAB
            for batchFile in ~{sep=" " batch_file_list}; do
                batchID=$( basename $batchFile ".txt")
                echo "${batchID}"
                
                bcftools view -S ${batchFile} --force-samples \
                ~{file_prefix}.annot.exon.0_01.imp.bcf -Ov -o \
                ~{file_prefix}.annot.exon.0_01.imp.${batchID}.vcf
                
                python /opt/vep/custom/vcf2tab.v2.py \
                ~{file_prefix}.annot.exon.0_01.imp.${batchID}.vcf -o \
                ~{file_prefix}.annot.exon.0_01.imp.${batchID}.tab

                rm ~{file_prefix}.annot.exon.0_01.imp.${batchID}.vcf*

                bgzip -c ~{file_prefix}.annot.exon.0_01.imp.${batchID}.tab > \
                ~{file_prefix}.annot.exon.0_01.imp.${batchID}.tab.bg.gz

                rm ~{file_prefix}.annot.0_01.imp.${batchID}.tab
            done
            
            rm ~{file_prefix}.annot.exon.0_01.imp.bcf*
        >>>
        output {
            Array[File] annot_bcf = glob("output/*.annot.bcf*")
            Array[File] tab_gz_file = glob("output/*.tab.bg.gz")
            #Array[File] g_vcf_index = glob("output/*.g.vcf.gz.tbi")
        }
         runtime {
            docker : vep_docker
            dx_timeout: "48H"
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
            }         }
}
    
