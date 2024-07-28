#------------------
# Description: Script to create EDU-GWAS-SNP matrix
#

require('data.table')
require('Matrix')
source('ukbb/regression.R')

getSnpDF <- function(){

    OBK_PATH = '/hps/nobackup/dunham/ai-uk-can/obkay_2022'

    edu_gene_snp_map = paste(OBK_PATH,'obk_gene_var_map_sort.txt',sep="/")
    edu_gwas_samples = paste(OBK_PATH,'bgen_comb/3925_snps_all_samples.txt',sep="/")

    gene_snp_data = fread(edu_gene_snp_map,sep='\t',stringsAsFactors=F,header=F)
    gwas_samples = fread(edu_gwas_samples,header=F,sep='_',stringsAsFactors=F)

    gene_snp_data = data.frame(gene_snp_data)
    gwas_samples_df = data.frame(gwas_samples)

    edu_gt_df = data.frame(gwas_samples_df$V1)
    colnames(edu_gt_df) = "f.eid_samples"
    edu_gt_df$f.eid_samples = as.character(edu_gt_df$f.eid_samples)

#----------------------
# Read gene-snp-map data
#-----------------------

    geneNumList = c()

    for (i in 1:dim(gene_snp_data)[1]){
    #for (i in 1:10){
        message(i)
        var_list = strsplit(gene_snp_data[i,2],',')[[1]]
        for (var_id in var_list){
            var_id_file = paste0(OBK_PATH,'/mat/',var_id,'.txt',collapse="")
            varData = fread(var_id_file,sep='\t',header=T,stringsAsFactors=F,check.names=F)
            varData[varData=="0/0"] = 0
            varData[varData=="0/1"] = 1
            varData[varData=="1/1"] = 2
            edu_gt_df = cbind(edu_gt_df,data.frame(varData))
        }
        geneNumList = c(geneNumList,gene_snp_data[i,3]) 
        #stop()

    }

    edu_gt_df_bkp = edu_gt_df
    edu_gt_df[edu_gt_df == "0/0"] = 0
    edu_gt_df[edu_gt_df == "0/1"] = 1
    edu_gt_df[edu_gt_df == "1/1"] = 2
    edu_gt_df[edu_gt_df == "./."] = 0

    outList = list()
    outList[[1]] = edu_gt_df
    outList[[2]] = geneNumList
    save(outList,file=paste(OBK_PATH,'dl_data/out_list_raw',sep="/"))
    #save(edu_gt_df,file=paste(OBK_PATH,'dl_data/edu_gt_df_raw',sep="/"))
    #write.table(geneNumList,file=paste(OBK_PATH,'dl_data/edu_snps_size.txt',sep="/"),row.names=F,col.names=F,quote=F)
#write.table(edu_ft_df[,-1],file=paste(OBK_PATH,'dl_data/genotype.csv',sep="/"),row.names=F,col.names=F,quote=F)

}

mergeSnpPTV <- function(outList){
   
    require('dplyr')
    #------------------------------------------------------
    # Merge Common variant and Rare variants for NEDU years
    #------------------------------------------------------

    OBK_PATH = '/hps/nobackup/dunham/ai-uk-can/obkay_2022'

    edu_gt_df = outList[[1]]
    geneNumList = outList[[2]]

    load('/hps/nobackup/dunham/ai-uk-can/analysis/qualf/mat/cgt.qualf.mat.frame_impact')
    nedu_frame_gene = 'NF1'
    nf1.ind = which(colnames(cgt.rm.df)=='NF1')
    nf1.frame = cgt.rm.df[,c(1,nf1.ind,235:261)]
    
    #Extract NF1 gene
    ukb.reg.df = nf1.frame
    ukb.lm.fit =  getRegression(ukb.reg.df[,-c(1:2)],"nedu","glm","ALL")
    ukb.nedu.res = residuals(ukb.lm.fit)
    ukb.reg.df$nedu_resid[is.na(ukb.reg.df$nedu)==FALSE] = ukb.nedu.res
    ukb.reg.nna.df = subset(ukb.reg.df,!is.na(nedu))
    ukb.nedu.resid = data.frame(ukb.reg.nna.df$nedu_resid)

    #Extract NEDU CV & RV Reg Pheno 
    edu_cv = edu_gt_df[edu_gt_df$f.eid_samples %in% ukb.reg.nna.df$f.eid_samples,] 
    message('-- After extracting NEDU CV and RV Reg Pheno')

    #Extract CV and NF1 genes
    edu.cv.nf1 = dplyr::inner_join(data.frame(edu_cv[,1:2]),data.frame(ukb.reg.nna.df[,1:2]),by=join_by(f.eid_samples == f.eid_samples))
    edu.cv.rv = cbind(edu_cv,edu.cv.nf1[,3])
    write.table(edu.cv.rv[,-1],sep=",",file=paste(OBK_PATH,'dl_data/genotype.csv',sep="/"),row.names=F,col.names=F,quote=F)
    message('-- After extracting CV and NF1 genes')

    #Extract CV and NEDU residualized 
    edu.cv.nedu_resid = dplyr::inner_join(data.frame(edu_cv[,1:2]),data.frame(ukb.reg.nna.df[,c(1,30)]),by=join_by(f.eid_samples == f.eid_samples))
    write.table(edu.cv.nedu_resid$nedu_resid,file=paste(OBK_PATH,'dl_data/phenotype.csv',sep="/"),row.names=F,col.names=F,quote=F)
    message('-- After extracting CV and NEDU residualized')
    
    geneNumList = c(geneNumList$V1,1)
    write.table(geneNumList,file=paste(OBK_PATH,'dl_data/snp_size.csv',sep="/"),row.names=F,col.names=F,quote=F)
    message('-- Gene Num List')


    cv_rv_nedu_list = list()
    cv_rv_nedu_list[[1]] = edu.cv.rv
    cv_rv_nedu_list[[2]] = edu.cv.nedu_resid
    cv_rv_nedu_list[[3]] = geneNumList
    save(cv_rv_nedu_list,file=paste(OBK_PATH,'dl_data/cv_rvframe_nedu_list',sep="/"))
    message('-- Aftersaving CV and RV NEDU objects')

}

quantNorm <- function(.data){
    data_sort <- apply(.data, 2, sort)
    row_means <- rowMeans(data_sort)
    data_sort <- matrix(row_means,nrow = nrow(data_sort), ncol = ncol(data_sort), byrow = TRUE)
    index_rank <- apply(.data, 2, order)
    normalized_data <- matrix(nrow = nrow(.data), ncol = ncol(.data))
    for(i in 1:ncol(.data)){
          normalized_data[,i] <- data_sort[index_rank[,i], i]
    }
    return(normalized_data)
}

mergeSnpDMV <- function(outList,cgt_type){
   
    require('dplyr')
    #------------------------------------------------------
    # Merge Common variant and Rare variants for NEDU years
    #------------------------------------------------------

    OBK_PATH = '/hps/nobackup/dunham/ai-uk-can/obkay_2022'

    edu_gt_df = outList[[1]]
    geneNumList = outList[[2]]

    if (cgt_type == "nedu") {
        load('/hps/nobackup/dunham/ai-uk-can/analysis/qualf/mat/cgt.qualf.mat.miss_impact')
        geneList = c("DPP6","ACTL6A")
        cov.ind = c(384:410)
    }else if (cgt_type == "fit"){
        geneList = c("RAI1","HECW2","SRP54","FBXO31","GNAS","HCCS","CDH2","NCKAP1")
    }else if(cgt_type == "tdi"){
        geneList = c("SMARCA5","DHX30","GLI2","PUM1","ATP6V1A","AGO2","ATRX","YWHAG","MAGEL2","PPP3CA","COL4A1","BCL11B","CHD2")
    }
        
    cgt.ind = which(colnames(cgt.rm.df)==cgt_type)

    #nedu_frame_gene = 'NF1'
    #nf1.ind = which(colnames(cgt.rm.df)=='NF1')
    gene.ind = which(colnames(cgt.rm.df) %in% geneList == TRUE)
    gene.dmv = cgt.rm.df[,c(1,gene.ind,384:cgt.ind)]
    gene.dmv.cov = cgt.rm.df[,cov.ind]
    
    #Residuals of regression covariates
    ukb.reg.df = gene.dmv
    ukb.lm.fit =  getRegression(gene.dmv.cov,"nedu","glm","ALL")
    ukb.nedu.res = residuals(ukb.lm.fit)
    ukb.reg.df$nedu_resid[is.na(ukb.reg.df$nedu)==FALSE] = ukb.nedu.res
    ukb.reg.nna.df = subset(ukb.reg.df,!is.na(nedu))
    ukb.nedu.resid = data.frame(ukb.reg.nna.df$nedu_resid)
    geneList.ind = which(colnames(ukb.reg.nna.df) %in% geneList == TRUE)
    reg.cgt.ind = which(colnames(ukb.reg.nna.df) %in% cgt_type ==TRUE)

    #Extract NEDU CV & RV Reg Pheno 
    edu_cv = edu_gt_df[edu_gt_df$f.eid_samples %in% ukb.reg.nna.df$f.eid_samples,] 
    message('-- After extracting NEDU CV and RV Reg Pheno')

    #Extract CV and EDU genes
    edu.cv.genes = dplyr::inner_join(data.frame(edu_cv[,1:2]),data.frame(ukb.reg.nna.df[,c(1,geneList.ind)]),by=join_by(f.eid_samples == f.eid_samples))

    edu.cv.rv = cbind(edu_cv,edu.cv.genes[,-c(1:2)])

    #Load EDU Top SNPs and Mappings
    load('/hps/nobackup/dunham/ai-uk-can/obkay_2022/ptv/dl_data/nedu/neduList')
    top.snp.gene.vars = neduList[[1]]
    top.gene.gt.mat = getGenotypeMat(edu.cv.rv,top.snp.gene.vars,geneList)
    geneNumList = top.snp.gene.vars$count
    geneNumList = c(geneNumList, rep(1,length(geneList)))

    message('-- After Loading NEDU List')

    write.table(top.gene.gt.mat[,-1],sep=",",file=paste(OBK_PATH,'dmv',cgt_type,'dl_data/genotype.csv',sep="/"),row.names=F,col.names=F,quote=F)
    message('-- After extracting CV and EDU genes')

    #Extract CV and NEDU residualized 
    edu.cv.nedu_resid = dplyr::inner_join(data.frame(top.gene.gt.mat[,1:2]),data.frame(ukb.reg.nna.df[,c(1,reg.cgt.ind,reg.cgt.ind+1)]),by=join_by(f.eid_samples == f.eid_samples))
    write.table(edu.cv.nedu_resid$nedu_resid,file=paste(OBK_PATH,'dmv',cgt_type,'dl_data/phenotype.csv',sep="/"),row.names=F,col.names=F,quote=F)
    message('-- After extracting CV and NEDU residualized')
    
    #geneNumList = c(geneNumList$V1,1)
    write.table(geneNumList,file=paste(OBK_PATH,'dmv',cgt_type,'dl_data/snp_size.csv',sep="/"),row.names=F,col.names=F,quote=F)
    message('-- Gene Num List')


    cv_rv_nedu_list = list()
    cv_rv_nedu_list[[1]] = top.gene.gt.mat
    cv_rv_nedu_list[[2]] = edu.cv.nedu_resid
    cv_rv_nedu_list[[3]] = geneNumList
    save(cv_rv_nedu_list,file=paste(OBK_PATH,'dmv',cgt_type,'dl_data/cv_rvframe_nedu_list',sep="/"))
    message('-- Aftersaving CV and RV NEDU objects')

    return(cv_rv_nedu_list)
}

getGenotypeMat <- function(edu.cv.rv, top.snp.gene.vars,geneList) {

    message('Createing Genotype Matrix for top SNPs')
    #Subroutine to create genotype matri w.r.t genes
    gt.df = data.frame(edu.cv.rv[,1])
    gt_col_names = c("f.eid_samples")

    for (i in 1:dim(top.snp.gene.vars)[1]){
        vars_strs = strsplit(top.snp.gene.vars$vars[i],",")[[1]]
        for(vars in vars_strs){
            ind = which(grepl(vars,colnames(edu.cv.rv))==TRUE)[1]
            gt.df = cbind(gt.df,as.numeric(edu.cv.rv[,ind]))
            gt_col_names = c(gt_col_names,vars)
        }
    }
    gene.ind = which(colnames(edu.cv.rv) %in% geneList==TRUE)
    gt.df = cbind(gt.df,edu.cv.rv[,gene.ind])
    gt_col_names = c(gt_col_names,colnames(edu.cv.rv)[gene.ind])
    colnames(gt.df) = gt_col_names

    return(gt.df)

}
