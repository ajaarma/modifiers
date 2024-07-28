###########################################################################
# Description: Script to merge hpc-200k and dnx-272k. Merging is done by
#   (a) Finding unique vars in HPC
#   (b) Finding unique vars in DNX
#   (c) Finding common vars in HPC & DNX
#   (d) Transcript Prioritzation (get the index)
#   (e) Combine the samples
#   (f) create three cases
#
# Author: Ajay A. Kumar (aak@ebi.ac.uk)
#
# Date: 26.12.2023
#
############################################################################

#args=(commandArgs(TRUE))

require('data.table')
require('R.utils')
source("/homes/aak/scripts/AI-UK/src/dv/rlib/TransPrior.R")
#sourceDirectory(paste("/homes/aak/scripts/AI-UK/src/dv/rlib",sep=""))
#gene_file = '/hps/nobackup/dunham/ai-uk-can/ukb200k/query_ID_cadd/GALC_ukb23156_c14_b20_v1_vars_uniq_samples.txt'

########## FUNCTIONS #############

checkHpcDnxFiles <- function(hpc_dir, dnx_dir){
    # Subroutine to compare difference of files in both directories
    
    hpcGeneList = list.files(hpc_dir)
    dnxGeneList = list.files(dnx_dir)

    geneFileDnxList = c()
    geneFileHpcList = c()

    for (ele in hpcGeneList){
        #strs = strsplit(ele,"_")[[1]]
        #gene_id = strs[1]
        #chrBatch = paste(strs[1],strs[3],strs[4],strs[5],sep="_")
        geneFileHpcList = c(geneFileHpcList,ele)
    }

    for (ele in dnxGeneList){
        #strs = strsplit(ele,"_")[[1]]
        #chrBatch = paste(strs[1],strs[3],strs[4],strs[5],sep="_")
        geneFileDnxList = c(geneFileDnxList,ele)
    }

    outList = list()
    outList[[1]] = geneFileDnxList
    outList[[2]] = geneFileHpcList

    return(outList)

}


processDnxOrHpc <- function(geneDnxFile){
    #Process only DNX-diff Files
    dnxVar = fread(geneDnxFile,sep='\t',header=T,stringsAsFactors=F,quote="")
    dnx_diff_df = data.frame(dnxVar)
    hpc_samples_list = rep('NA',dim(dnx_diff_df)[1])
    hpc_samples_count_list = rep(0,dim(dnx_diff_df)[1])
    SumCount = as.vector(dnx_diff_df[,'Count'])
    dnx_diff_df = cbind(dnx_diff_df,hpc_samples_list,hpc_samples_count_list,SumCount)

    return(dnx_diff_df)
}

processDnxHpc <- function(geneDnxFile,geneHpcFile){
    #Process DNX and HPC Files

    dnxVar = fread(geneDnxFile,sep='\t',header=T,stringsAsFactors=F,quote="")
    hpcVar = fread(geneHpcFile,sep='\t',header=T,stringsAsFactors=F,quote="")

    dnx_hpc_vars = intersect(dnxVar$VAR_KEY,hpcVar$VAR_KEY)
    dnx_diff_vars = setdiff(dnxVar$VAR_KEY,hpcVar$VAR_KEY)
    hpc_diff_vars = setdiff(hpcVar$VAR_KEY,dnxVar$VAR_KEY)

    dnx_hpc_df = subset(dnxVar,VAR_KEY %in% dnx_hpc_vars)
    dnx_diff_df = subset(dnxVar,VAR_KEY %in% dnx_diff_vars)
    hpc_diff_df = subset(hpcVar,VAR_KEY %in% hpc_diff_vars)
    #dnx_diff_df

    hpc_samples_list = c()
    hpc_samples_count_list = c()
    dnx_diff_samples = c()
    hpc_diff_samples = c()
    SumCount = c()

    #Processing dnx-hpc-df
    dnx_hpc_df = data.frame(dnx_hpc_df)
    for (i in 1:dim(dnx_hpc_df)[1]){
        var_key = dnx_hpc_df[i,'VAR_KEY']
        dnx_samples = dnx_hpc_df[i,'UniqSamples']
        dnx_samples_count = dnx_hpc_df[i,'Count']

        hpc_samples = as.vector(subset(hpcVar,VAR_KEY==var_key)[,'UniqSamples']$UniqSamples)
        hpc_samples_count = as.vector(subset(hpcVar,VAR_KEY==var_key)[,'Count']$Count)
        hpc_samples_count_list = c(hpc_samples_count_list,hpc_samples_count)
        hpc_samples_list = c(hpc_samples_list,hpc_samples)
        SumCount = c(SumCount,dnx_samples_count+hpc_samples_count)
    }
    dnx_hpc_df = cbind(dnx_hpc_df,data.frame(hpc_samples_list),hpc_samples_count_list,SumCount)

    #Processing dnx-diff-hpc
    dnx_diff_df = data.frame(dnx_diff_df)
    hpc_samples_list = rep('NA',dim(dnx_diff_df)[1])
    hpc_samples_count_list = rep(0,dim(dnx_diff_df)[1])
    SumCount = as.vector(dnx_diff_df[,'Count'])
    dnx_diff_df = cbind(dnx_diff_df,hpc_samples_list,hpc_samples_count_list,SumCount)

    #Processing hpc-diff-df
    hpc_diff_df = data.frame(hpc_diff_df)
    hpc_samples_list = hpc_diff_df[,'UniqSamples']
    hpc_samples_count_list = hpc_diff_df[,'Count']
    SumCount = as.vector(hpc_diff_df[,'Count'])
    hpc_diff_df = cbind(hpc_diff_df,hpc_samples_list,hpc_samples_count_list,SumCount)

    #dnx_hpc_merge_df = rbind(dnx_hpc_df,dnx_diff_df,hpc_diff_df,fill=T)
    tryCatch(
             {
                dnx_hpc_merge_df = rbind(dnx_hpc_df,dnx_diff_df,hpc_diff_df)
                outList = list()
                outList[[1]] = dnx_hpc_vars
                outList[[2]] = dnx_diff_vars
                outList[[3]] = hpc_diff_vars
                outList[[4]] = dnxVar
                outList[[5]] = hpcVar
                outList[[6]] = dnx_hpc_merge_df 
             
            },error = function(e){
                print(e)
             },finally = {
                
                dnx_hpc_merge_df = rbindlist(list(dnx_hpc_df,dnx_diff_df,hpc_diff_df),fill=T)
                outList = list()
                outList[[1]] = dnx_hpc_vars
                outList[[2]] = dnx_diff_vars
                outList[[3]] = hpc_diff_vars
                outList[[4]] = dnxVar
                outList[[5]] = hpcVar
                outList[[6]] = dnx_hpc_merge_df 
             }
    )
    return(outList)
}

mergeGeneFiles <- function(query_id_dir,out_dir){
    
    #Subroutine to merge gene files in DNX/HPC directories
    geneFileList = list.files(query_id_dir)
    geneList = c()

    for (ele in geneFileList){
        strs = strsplit(ele,"_")[[1]]
        gene_id = strs[1]
        geneList = c(geneList,gene_id)
    }
   
    geneList = unique(geneList)

    for (gene_id in geneList){
        geneFileList.2 = list.files(query_id_dir,pattern=paste0("^",gene_id,"_",collapse=""))
        outList = list()
        print(geneFileList.2)

        for (i in 1:length(geneFileList.2)){
            geneFile = paste(query_id_dir,geneFileList.2[i],sep='/')
            strs = strsplit(geneFile,'_')[[1]]
            rawVar = fread(geneFile,sep='\t',header=T,stringsAsFactors=F,quote="")
            outList[[i]] = rawVar
        }
        strs = strsplit(geneFileList.2[1],'_')[[1]]
        outGeneFile = paste(out_dir,paste(strs[1],'vars_merge_uniq_samples.txt',sep="_"),sep='/')
        #merge.df = do.call("rbind",outList,fill=TRUE)
        merge.df = rbindlist(outList,fill=T)
        write.table(unique(merge.df),sep='\t',file=outGeneFile,row.names=F,col.names=T,quote=F)
        #return(merge.df)
    }
}

trans_map = '/nfs/research/dunham/resources/ensembl/grch38/ensBioMart_grch38_v100_ENST_lengths_200510.txt'
transData = fread(trans_map,sep='\t',header=T,stringsAsFactors=F,quote="")
transData = data.frame(transData)
colnames(transData) = c("Gene_StableID","Transcript_StableID","Transcript_length","Transcript_type","Gene_name","Gene_type")

hpc_cadd_dir = '/hps/nobackup/dunham/ai-uk-can/ukb200k/query_ID_cadd'
dnx_cadd_dir = '/hps/nobackup/dunham/ai-uk-can/dnx_280k/query_ID_cadd'

hpc_merge_dir = '/hps/nobackup/dunham/ai-uk-can/ukb200k/query_ID_merge'
dnx_merge_dir = '/hps/nobackup/dunham/ai-uk-can/dnx_280k/query_ID_merge'

outSetDir = '/hps/nobackup/dunham/ai-uk-can/analysis/hpc_dnx_ukbb'
outPriorDir = '/hps/nobackup/dunham/ai-uk-can/analysis/hpc_dnx_prior_ukbb'
######### Merge files per Gene : DNX/HPC batches ###########

#mergeGeneFiles(hpc_cadd_dir,hpc_merge_dir) #Comment out to produce the output in merge directory
#mergeGeneFiles(dnx_cadd_dir,dnx_merge_dir)

#stop()
########## SET analysis of DNX-HPC files #################

outList = checkHpcDnxFiles(hpc_merge_dir,dnx_merge_dir)
geneFileDnxList = outList[[1]]
geneFileHpcList = outList[[2]]
dnx_hpc = intersect(geneFileDnxList,geneFileHpcList)
dnx_diff = setdiff(geneFileDnxList,geneFileHpcList)
hpc_diff = setdiff(geneFileHpcList,geneFileDnxList)
    
logFile = paste(outPriorDir,'.logFile.log',sep='/')

sink(logFile)
########## Processing DNX-DIFF-PRIOR ####################
for (ele in dnx_diff){
    gene_id = strsplit(ele,'_')[[1]][1]
    message('Processing DNX-DIFF-Prior: ',gene_id)
    cat('Processing DNX-DIFF-Prior: ',gene_id,'\n')
    geneDnxFile = paste(dnx_merge_dir,ele,sep="/")
    outList.merge = processDnxOrHpc(geneDnxFile)
    outFilePath = paste(outSetDir,ele,sep='/')
    write.table(outList.merge,sep="\t",file=outFilePath,row.names=F,col.names=T,quote=F)
    
    varData.prior = getTransPrior(outList.merge,transData,gene_id)
    outPriorFile = paste(outPriorDir,paste(gene_id,'vars_prior.txt',sep="_"),sep='/')
    write.table(varData.prior,sep='\t',file=outPriorFile,row.names=F,col.names=T,quote=F)
}

########## Processing HPC-DIFF-PRIOR ####################

for (ele in hpc_diff){
    gene_id = strsplit(ele,'_')[[1]][1]
    message('Processing HPC-DIFF-Prior: ',gene_id)
    cat('Processing HPC-DIFF-Prior: ',gene_id,'\n')
    
    geneDnxFile = paste(hpc_merge_dir,ele,sep="/")
    outList.merge = processDnxOrHpc(geneDnxFile)
    outFilePath = paste(outSetDir,ele,sep='/')
    write.table(outList.merge,sep="\t",file=outFilePath,row.names=F,col.names=T,quote=F)

    varData.prior = getTransPrior(outList.merge,transData,gene_id)
    outPriorFile = paste(outPriorDir,paste(gene_id,'vars_prior.txt',sep="_"),sep='/')
    write.table(varData.prior,sep='\t',file=outPriorFile,row.names=F,col.names=T,quote=F)
}

########## Processing DNX-HPC-PRIOR ####################

for (ele in dnx_hpc){
    gene_id = strsplit(ele,'_')[[1]][1]
    message('Processing DNX-HPC-Prior: ',gene_id)
    cat('Processing DNX-HPC-Prior: ',gene_id,'\n')
    
    geneDnxFile = paste(dnx_merge_dir,ele,sep="/")
    geneHpcFile = paste(hpc_merge_dir,ele,sep="/")
    outList.merge = processDnxHpc(geneDnxFile,geneHpcFile)
    outFilePath = paste(outSetDir,ele,sep='/')
    write.table(outList.merge[[6]],sep="\t",file=outFilePath,row.names=F,col.names=T,quote=F)

    varData.prior = getTransPrior(outList.merge[[6]],transData,gene_id)
    outPriorFile = paste(outPriorDir,paste(gene_id,'vars_prior.txt',sep="_"),sep='/')
    write.table(varData.prior,sep='\t',file=outPriorFile,row.names=F,col.names=T,quote=F)
 
}

sink()
