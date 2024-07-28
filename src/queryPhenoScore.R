require('Matrix')
require('data.table')
require('R.utils')
require('dplyr')
require('R.utils')

source('ukbb/prs.R')
prsFile = '/hps/nobackup/dunham/ai-uk-can/obkay_2022/plink_conv/maf_01/PRS.sscore'
res_05 = '/hps/nobackup/dunham/ai-uk-can/analysis/results/fdr/out_vars_0_05.txt'

resDat = read.table(res_05,sep='\t',header=T)
prs.df.slice = getPRScut(prsFile)

cgtList = c('qualf','fit','tdi')

for (cgt_type in cgtList){
    if(cgt_type == 'qualf'){

        message('Processing NEDU')
        #PTVs
        load('/hps/nobackup/dunham/ai-uk-can/analysis/qualf/mat/cgt.qualf.mat.frame_impact')
        cgt.gene.df = subset(resDat,cgt_imp_class == cgt_type & var_imp_class=='frame')
        cgt.gene.neg = subset(cgt.gene.df,Lasso.Est <0)$Predictors
        gene.ind = which(colnames(cgt.rm.df)==cgt.gene.neg)
        fid.ind = which(colnames(cgt.rm.df)=='f.eid_samples')
        nedu.ind = which(colnames(cgt.rm.df)=='nedu')
        cgt.gene.nedu = cgt.rm.df[,c(fid.ind,nedu.ind,gene.ind)]
        cgt.gene.nedu.ca = cgt.gene.nedu[cgt.gene.nedu[,3]==1,]
        cgt.gene.nedu.ca.sd1 = subset(cgt.gene.nedu.ca,nedu > 1 | nedu < -1)

        ### Merge with PRS

        nedu.sd1.prs = dplyr::left_join(cgt.gene.nedu.ca.sd1,prs.df.slice,by = join_by(f.eid_samples == FID))
        nedu.prs = dplyr::left_join(cgt.gene.nedu.ca,prs.df.slice,by = join_by(f.eid_samples == FID))


    }else if(cgt_type == 'fit'){
        message('Processing FIT')

        load('/hps/nobackup/dunham/ai-uk-can/analysis/fit/mat/cgt.fit.mat.frame_impact')
        cgt.gene.df = subset(resDat,cgt_imp_class == cgt_type & var_imp_class=='frame')
        cgt.gene.neg = subset(cgt.gene.df,Lasso.Est <0)$Predictors
        gene.ind = which(colnames(cgt.rm.df)==cgt.gene.neg)
        fid.ind = which(colnames(cgt.rm.df)=='f.eid_samples')
        fit.ind = which(colnames(cgt.rm.df)=='fit')
        cgt.gene.fit = cgt.rm.df[,c(fid.ind,fit.ind,gene.ind)]
        cgt.gene.fit.ca = cgt.gene.fit[cgt.gene.fit[,3]==1,]
        cgt.gene.fit.ca.sd1 = subset(cgt.gene.fit.ca,fit > 1 | fit < -1)

        ### Merge with PRS

        fit.sd1.prs = dplyr::left_join(cgt.gene.fit.ca.sd1,prs.df.slice,by = join_by(f.eid_samples == FID))
        fit.prs = dplyr::left_join(cgt.gene.fit.ca,prs.df.slice,by = join_by(f.eid_samples == FID))
    
    }else if (cgt_type == 'tdi'){
        message('Processing TDI')
        load('/hps/nobackup/dunham/ai-uk-can/analysis/tdi/mat/cgt.tdi.mat.frame_impact')
        cgt.gene.df = subset(resDat,cgt_imp_class == cgt_type & var_imp_class=='frame' & filt_imp_class=='impact')
        cgt.gene.neg = subset(cgt.gene.df,Lasso.Est >0)$Predictors
        gene.ind = which(colnames(cgt.rm.df)==cgt.gene.neg[1])
        fid.ind = which(colnames(cgt.rm.df)=='f.eid_samples')
        tdi.ind = which(colnames(cgt.rm.df)=='tdi')
        cgt.gene.tdi = cgt.rm.df[,c(fid.ind,tdi.ind,gene.ind)]
        cgt.gene.tdi.ca = cgt.gene.tdi[cgt.gene.tdi[,3]==1,]
        cgt.gene.tdi.ca.sd1 = subset(cgt.gene.tdi.ca,tdi > 1 | tdi < -1)

        ### Merge with PRS

        tdi.sd1.prs = dplyr::left_join(cgt.gene.tdi.ca.sd1,prs.df.slice,by = join_by(f.eid_samples == FID))
        tdi.prs = dplyr::left_join(cgt.gene.tdi.ca,prs.df.slice,by = join_by(f.eid_samples == FID))
    
    }

}


df = data.frame(matrix(ncol=7,nrow=5))
  colnames(df) = c('qntl','NF1_NEDU_ALL','NF1_NEDU_SD1','MAPK1_FIT_ALL','MAPK1_FIT_SD1','SPTBN1_TDI_ALL','SPTBN1_TDI_SD1')
  df$qntl = c('0-10%',"20-40%","40-60%","60-80%","80-100%")
  
  #### NEDU ######
  df$NF1_NEDU_ALL = c(18,17,21,21,14)
  df$NF1_NEDU_SD1 = c(9,7,13,11,4)  
  
  #### FIT #####
  df$MAPK1_FIT_ALL = c(2,2,0,2,1)
  df$MAPK1_FIT_SD1 = c(0,1,0,1,0)
  
  #### TDI #####
  df$SPTBN1_TDI_ALL = c(1,2,5,1,3)
  df$SPTBN1_TDI_SD1 = c(0,0,4,0,1)
  
  df.melt = melt(df,"qntl")
  
  p = ggplot(df.melt,aes(x=qntl,y=value,fill=variable))+
    geom_bar(stat="identity",position="dodge")+
    #ggtitle("NEDU individuals per quintile")+
    xlab("PRS Qantiles")+ylab("No. of Samples")+
    theme(axis.title.x = element_text(size = 14,vjust=-1.0,face="bold"), 
          axis.title.y = element_text(size = 14,hjust=0.5,face="bold"),
          text=element_text(size=12,face="bold"),
          plot.title=element_text(size=16,hjust=0.5,face="bold"))+
    theme(legend.title=element_blank())
  #print(p)
  #ggsave(p,file="/Users/ajaykumar/data/AI-CAN-UK/ukbb/results/prs/gene_samples_cgt_ptv.png")


