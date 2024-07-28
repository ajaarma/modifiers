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
task mafImpTab {
    input {
        File anno_bcf
        File anno_bcf_index
        File norm_bcf
        File norm_bcf_index
        File maf_imp_bcf
        Array[File] batch_file_list
        String file_prefix = basename(norm_bcf,".bcf")
        String imp_prefix = basename(maf_imp_bcf,".bcf")
        String vep_docker
    }    
    
    command <<<
        
        export PATH=/opt/vep/bcftools-1.9/bin:$PATH
        export BCFTOOLS_PLUGINS=/opt/vep/bcftools-1.9/plugins
        export PATH=/opt/vep/custom:$PATH

        #generate maf Impact index
        bcftools index ~{maf_imp_bcf}
        #mv ~{file_prefix}.csi ~/inputs/

        #Merge wGT-annot & AC0.norm
        echo "Merging wGT and ACo.norm"
        bcftools annotate -a ~{anno_bcf} \
        ~{norm_bcf} -c ID,FILTER,INFO \
        -Ob -o ~{file_prefix}.annot.merge.bcf
        bcftools index ~{file_prefix}.annot.merge.bcf

        rm ~{norm_bcf}
        rm ~{norm_bcf_index}

        #Intersect betwen maf-imp & norm.annot file
        echo "Intersecting bcf files"
        bcftools isec -n=2 ~{maf_imp_bcf} \
        ~{file_prefix}.annot.merge.bcf -p ./ -Ob

        #Rename isec file 
        mv 0001.bcf ~{imp_prefix}.merge.bcf
        mv 0001.bcf.csi ~{imp_prefix}.merge.bcf.csi

        #Convert VCF 2 TAB
        for batchFile in ~{sep=" " batch_file_list}; do
            batchID=$( basename $batchFile ".txt")
            echo "${batchID}"
            bcftools view -S ${batchFile} --force-samples \
            ~{imp_prefix}.merge.bcf -Ov -o \
            ~{imp_prefix}.merge.${batchID}.vcf

            python /opt/vep/custom/vcf2tab.v2.py \
            ~{imp_prefix}.merge.${batchID}.vcf -o \
            ~{imp_prefix}.merge.${batchID}.tab

            rm ~{imp_prefix}.merge.${batchID}.vcf*

            bgzip -c ~{imp_prefix}.merge.${batchID}.tab > \
            ~{imp_prefix}.merge.${batchID}.tab.bg.gz
            rm ~{imp_prefix}.merge.${batchID}.tab
        done
    >>>
    
    output{
        File anno_merge_bcf = "~{file_prefix}.annot.merge.bcf"
        File anno_merge_index = "~{file_prefix}.annot.merge.bcf.csi"
        File imp_merge_bcf = "~{imp_prefix}.merge.bcf"
        File imp_merge_index = "~{imp_prefix}.merge.bcf.csi"
        File tab_1 = "~{imp_prefix}.merge.batch1.tab.bg.gz"
        File tab_2 = "~{imp_prefix}.merge.batch2.tab.bg.gz"
        File tab_3 = "~{imp_prefix}.merge.batch3.tab.bg.gz"
        File tab_4 = "~{imp_prefix}.merge.batch4.tab.bg.gz"
    }
    runtime {
        docker : vep_docker
        dx_timeout: "48H"
        #dx_instance_type : "mem2_ssd1_v2_x2"
    }
}

