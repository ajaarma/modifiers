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

task sortByExon {
    input {
        File exon_region = "dx://500K-WES:/resources/ensembl/ensembl_hg38.default.exons.50.sorted.merged.bed.gz"
        File annot_bcf
        File annot_bcf_index

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
        docker : "dx://project-GV63Xq0Jkqky3PqPj7Jp1VPZ:/resources/soft/bcftools-1.19.tar"
        dx_timeout : "48H"
        dx_instance_type : "mem3_ssd2_v2_x2"
    }
}

