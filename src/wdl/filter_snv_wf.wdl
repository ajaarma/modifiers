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
task annoMergeExon {
        
        input {
            #Input parameters from command line
            File input_snv_file
            File input_snv_index
            File anno_gt_file
            File anno_gt_index
             
            # Running parameters
            String bcf_docker

            String file_prefix = basename(input_snv_file,".bcf")
        }

        command <<<

            #Commented out because pipe command giving error

            #mkdir -p output/
            #cd output/
          
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
            File annot_merge = "~{file_prefix}.annot.bcf"
            File annot_merge_index = "~{file_prefix}.annot.bcf.csi"

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
            }         }
}

task sortByExon {
    input {
        File exon_region
        File annot_bcf
        File annot_bcf_index

        String bcf_docker
        String file_prefix = basename(annot_bcf,".bcf")
    }    
    
    command <<<
        
        #Extract exonic region
        echo "Extract exonic regions"
        bcftools view -R ~{exon_region} \
        ~{annot_bcf} -Ob -o \
        ~{file_prefix}.exon.bcf


        #Sort by exonic region coordinates
        bcftools view ~{file_prefix}.exon.bcf | \
        bcftools sort -m 1G -T /home/dnanexus/output | \
        bcftools norm -D | bcftools view -Ob -o \
        ~{file_prefix}.exon.sorted.bcf

        mv ~{file_prefix}.exon.sorted.bcf \ 
        ~{file_prefix}.exon.bcf
        bcftools index ~{file_prefix}.exon.bcf

    >>>
    
    output {
        File exon_bcf = "~{file_prefix}.exon.bcf"
        File exon_index = "~{file_prefix}.exon.bcf.csi"
    }
    runtime {
        docker : bcf_docker
        dx_timeout : "48H"
        dx_instance_type : "mem3_ssd2_v2_x2"
    }
}

task mafImpactTab {
    input {
        String exon_bcf
        String exon_bcf_index
        String vep_docker
        String chrNum
        Array[File] batch_file_list
        String file_prefix = basename(exon_bcf,".bcf")
    }    
    
    command <<<
        export PATH=/opt/vep/bcftools-1.9/bin:$PATH
        export BCFTOOLS_PLUGINS=/opt/vep/bcftools-1.9/plugins
        export PATH=/opt/vep/custom:$PATH

        if [ $chrNum != "chrY"]; then
            #Filter against gnomAD MAF <=0.01 or 1%
            echo "Filter against gnomAD MAF <=0.01"
            bcftools view -Ob -e 'GNOMADg4_AF>0.01' \
            ~{exon_bcf} | \
            bcftools view -Ob -e 'GNOMADe4_AF>0.01' > \
            ~{file_prefix}.0_01.bcf
        else
            echo "Filter against gnomAD MAF <=0.01 for chrY"
            bcftools view -Ob -e 'GNOMADe4_AF>0.01' \
            ~{exon_bcf} > ~{file_prefix}.0_01.bcf
        fi

        bcftools index ~{file_prefix}.0_01.bcf

        #Filter by Impact
        echo "Filter by Impact"
        bcftools view ~{file_prefix}.0_01.bcf | \
        grep -E
        '^#|^MT|MODERATE|HIGH|splice_donor_5th_base_variant|splice_donor_0th_base_variant|splice_donor_region_variant|splice_polypyrimidine_tract_variant|HGMD|CLNSIG' | \
        bcftools view -Ob -o ~{file_prefix}.0_01.imp.bcf
        bcftools index ~{file_prefix}.0_01.imp.bcf

        #Convert VCF 2 TAB
        for batchFile in ~{sep=" " batch_file_list}; do
            batchID=$( basename $batchFile ".txt")
            echo "${batchID}"
            bcftools view -S ${batchFile} --force-samples \
            ~{file_prefix}.0_01.imp.bcf -Ov -o \
            ~{file_prefix}.0_01.imp.${batchID}.vcf

            python /opt/vep/custom/vcf2tab.v2.py \
            ~{file_prefix}.0_01.imp.${batchID}.vcf -o \
            ~{file_prefix}.0_01.imp.${batchID}.tab

            rm ~{file_prefix}.0_01.imp.${batchID}.vcf*

            bgzip -c ~{file_prefix}.0_01.imp.${batchID}.tab > \
            ~{file_prefix}.0_01.imp.${batchID}.tab.bg.gz
            rm ~{file_prefix}.0_01.imp.${batchID}.tab
        done
    >>>
    
    output{
        File imp_bcf = "~{file_prefix}.0_01.imp.bcf"
        File imp_bcf_index = "~{file_prefix}.0_01.imp.bcf.csi"
        File tab_1 = "~{file_prefix}.0_01.imp.batch1.tab.bg.gz"
        File tab_2 = "~{file_prefix}.0_01.imp.batch2.tab.bg.gz"
        File tab_3 = "~{file_prefix}.0_01.imp.batch3.tab.bg.gz"
    }
    runtime {
        docker : vep_docker
        dx_timeout: "48H"
        dx_instance_type : "mem2_ssd1_v2_x2"
    }
}

workflow filterSnvWF {
        
        input{
            File input_snv_file
            File input_snv_index
            File anno_gt_file
            File anno_gt_index
            File exon_region
            String chrNum
            Array[File] batch_file_list

            #Default parameters
            String vep_docker
            String bcf_docker
        }
        String file_prefix = basename(input_snv_file,".bcf")

        call annoMergeExon {
            input: 
                input_snv_file = input_snv_file,
                input_snv_index = input_snv_index,
                anno_gt_file = anno_gt_file,
                anno_gt_index = anno_gt_index,
                bcf_docker = bcf_docker
        }
        call sortByExon {
            input:
                exon_region = exon_region,
                annot_bcf = annoMergeExon.annot_merge,
                annot_bcf_index = annoMergeExon.annot_merge_index,
                bcf_docker = bcf_docker
        }
        call mafImpactTab {
            input:
                exon_bcf = sortByExon.exon_bcf,
                exon_bcf_index = sortByExon.exon_index,
                batch_file_list = batch_file_list,
                chrNum = chrNum,
                vep_docker = vep_docker
        }
        output {
            File annot_merge = annoMergeExon.annot_merge
            File annot_merge_index = annoMergeExon.annot_merge_index
            File anno_merge_exon = sortByExon.exon_bcf
            File anno_merge_exon_index = sortByExon.exon_index
            File imp_bcf = mafImpactTab.imp_bcf
            File imp_bcf_index = mafImpactTab.imp_bcf_index
            File tab_1 = mafImpactTab.tab_1
            File tab_2 = mafImpactTab.tab_2
            File tab_3 = mafImpactTab.tab_3
        }
}
