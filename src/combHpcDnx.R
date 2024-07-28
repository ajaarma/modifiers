#########################################################################
# Description: Script to combine hpc-200k and dnx-272k filtered variants
#
# Email: aak@ebi.ac.uk
# Date: 25.12.2023
#
########################################################################


args=(commandArgs(TRUE))

require('data.table')
require('R.utils')


combFields <- function(df,anal_type){
    # Subroutine to get unique combination of the df field

    df$CHROM = gsub(" ","",df$CHROM)
    df$POS = gsub(" ","",df$POS)
    df$REF = gsub(" ","",df$REF)
    df$ALT = gsub(" ","",df$ALT)

    if (anal_type =='hpc') {
        gene.key = apply(data.frame(df),1,function(x){
                            strs = paste(x[2],x[3],x[5],x[6],sep="_")
                            return(strs)
                    })
    }else if (anal_type == 'cadd'){
        gene.key = apply(data.frame(df),1,function(x){
                            strs = paste(x[1],x[2],x[4],x[5],sep="_")
                            return(strs)
                    })
    
    
    }
    df.2 = cbind(gene.key,df)
    colnames(df.2) = c("VAR_KEY",colnames(df))

    return(df.2)
}

mapCadd <- function(ukb_dir,cadd_dir,anal_type){
    #Subroutine to map cadd score to query ID genes 

    query_id_dir = paste(ukb_dir,'query_ID_ukbb/uniq/',sep="/")
    outDir = paste(ukb_dir,'query_ID_cadd',sep='/')
    geneFileList = list.files(query_id_dir,pattern=".txt")

    for (ele in geneFileList){
        strs = strsplit(ele,"_")[[1]]
        chrBatch = paste(strs[2],strs[3],strs[4],strs[5],sep="_")
        chrNum = paste0('chr',strsplit(strs[3],"c")[[1]][2],collapse="")
        
        geneFilePath = paste0(query_id_dir,'/',ele,collapse="")
        gene.df = fread(geneFilePath,sep='\t',header=T,stringsAsFactors=F,quote="")

        if (anal_type=='hpc'){
            tabFile = 'tmpFile_tags.AC0.norm.anno.sites.merge.vep.tab'
        }else if (anal_type == 'dnx'){
            tabFile = 'tmpFile_tags.AC0.norm.anno.sites.exon.vep.tab'
        }

        caddFilePath = paste(cadd_dir,chrNum,chrBatch,tabFile,sep='/')

        cadd.df = fread(caddFilePath,sep='\t',header=T,stringsAsFactors=F,quote=F,fill=T)
       
        cat(ele,'\t',caddFilePath,'\n')
        outGeneFile = paste(outDir,ele,sep='/')
      
        gene.df.2 = combFields(gene.df,'hpc')
        cadd.df.2 = combFields(cadd.df,'cadd')

        #gene.cadd.tmp = cadd.df.2[cadd.df.2$VAR_KEY %in% gene.df.2$VAR_KEY,]
       
        anno_cadd_phred_list = c()
        anno_cadd_raw_list = c()

        for (i in 1:dim(gene.df.2)[1]){
            var_id = gene.df.2$VAR_KEY[i]
            anno_cadd = as.vector(subset(cadd.df.2,VAR_KEY==var_id)[,'ANNO_CADD_PHRED'])
            anno_raw = as.vector(subset(cadd.df.2,VAR_KEY==var_id)[,'ANNO_CADD_RAW'])
            #anno_cadd_phred_list = c(anno_cadd_phred_list,anno_cadd)
            #anno_cadd_raw_list = c(anno_cadd_raw_list,anno_raw)
            
            gene.df.2[i,'ANNO_CADD_PHRED'] = anno_cadd
            gene.df.2[i,'ANNO_CADD_RAW'] = anno_raw
        }

        #gene.df.2$ANNO_CADD_PHRED = anno_cadd_phred_list
        #gene.df.2$ANNO_CADD_RAW = anno_cadd_raw_list

        #return(gene.df.2)
        write.table(gene.df.2,sep='\t',file=outGeneFile,row.names=F,col.names=T,quote=F)
    }

}

######## HPC and DNX UKBB ############ 
hpc_ukb_dir = '/hps/nobackup/dunham/ai-uk-can/ukb200k' #args[[1]]
dnx_ukb_dir = '/hps/nobackup/dunham/ai-uk-can/dnx_280k' #args[[2]]

####### HPC & DNX query ID genes ##########
hpc_qid_dir = paste(hpc_ukb_dir,'query_ID_ukbb/uniq/',sep="")
dnx_qid_dir = paste(dnx_ukb_dir,'query_ID_ukbb/uniq/',sep="")

###### CADD Anno director ##############
hpc_cadd_dir = paste(hpc_ukb_dir,'cadd_anno_200k',sep="/")
dnx_cadd_dir = paste(dnx_ukb_dir,'dnx_cadd_b272k',sep="/")


######### Data processing #####

#outList = mapCadd(hpc_ukb_dir,hpc_cadd_dir,'hpc')
outList = mapCadd(dnx_ukb_dir,dnx_cadd_dir,'dnx')





