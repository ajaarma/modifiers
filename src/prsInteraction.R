require('data.table')
require('Hmisc')
require('moments')
require('dplyr')
require('R.utils')
sourceDirectory(paste("/homes/aak/scripts/AI-UK/src/ukbb/",sep=""))

#prsFile = '/hps/nobackup/dunham/ai-uk-can/impute/plink_conv/maf_005/PRS.sscore'
#prsFile = '/hps/nobackup/dunham/ai-uk-can/impute/plink_conv/maf_01/PRS.sscore'
#prsFile = '/hps/nobackup/dunham/ai-uk-can/obkay_2022/plink_conv/maf_01/PRS.sscore'
prsFile = '/hps/nobackup/dunham/ai-uk-can/obkay_2022/plink_conv/maf_005/PRS.sscore'
run_type = 'medu'
anal_type = 'frame'
reg_type = 'glm'

getPRScut <- function(prsFile) {
    prs.df = read.table(prsFile,sep='\t',header=F)
    prs.df$V5 = scale(prs.df$V5)
    quantiles = quantile(prs.df$V5,prob=seq(0,1,0.2))
    prs.df$slice = cut(prs.df$V5,breaks=c(min(prs.df$V5),-0.46335800,-0.06315392,0.26677819,0.64186147,max(prs.df$V5)),
                                    labels = c(1,2,3,4,5), include.lowest = TRUE)

    #prs.df.sort$slice = cut2(prs.df.sort$V5,cuts=as.numeric(quantiles),include.lowest=T,type=4)
    colnames(prs.df) = c('FID','IID','ALLELE_CT','NAMED_ALLELE_DOSAGE_SUM','PRS','SLICE')
    prs.df$FID = as.character(prs.df$FID)
    #prs.df.sort$PRS = scale(prs.df.sort$PRS)
    prs.df$PRS = prs.df$PRS

    return(prs.df)
}
getPtvMissPli <- function() {
    #----------------------------------------------------
    # Subroutine to get PTVs and Missense carrier status
    #----------------------------------------------------

    frame.sum.df = mergeFrameMissVars('frame')
    miss.sum.df = mergeFrameMissVars('miss')
    frame.miss.pli = dplyr::full_join(frame.sum.df,miss.sum.df,by='f.eid_samples')
    #frame.miss.pli = frame.sum.df
    colnames(frame.miss.pli)[c(2,3)] = c('gene.sum.frame','gene.sum.miss')
    frame.miss.pli$gene.sum.frame[frame.miss.pli$gene.sum.frame>0] = 1
    frame.miss.pli$gene.sum.frame[is.na(frame.miss.pli$gene.sum.frame)] = 0
    frame.miss.pli$gene.sum.miss[frame.miss.pli$gene.sum.miss>0] = 1
    frame.miss.pli$gene.sum.miss[is.na(frame.miss.pli$gene.sum.miss)] = 0

    return(frame.miss.pli)
}

glmnet_result <- function(resList){
    coef.lnz = resList$fitCV$coef.nz
    fm = resList[[3]]
    res.df = data.frame(summary(fm)$coefficients)
    res.df.2 = data.frame(cbind(rownames(res.df),res.df))
    colnames(res.df.2) = c("Predictors","GLM-Est","StdError","TValue","PValue")

    res.lnz.glm = cbind(coef.lnz,res.df.2[match(coef.lnz$Predictors,res.df.2$Predictors),])
    res.lnz.glm.2 = res.lnz.glm[,-3]
    colnames(res.lnz.glm.2) = c("Predictors","Lasso-Est","GLM-Est","GLM-StdError","TValue","PValue")
    res.lnz.glm.sort = res.lnz.glm.2[order(res.lnz.glm.2$PValue,decreasing=F),]

    return(res.lnz.glm.sort)

}
frame.miss.pli = getPtvMissPli()

#------------------------------------------
# Merge other covariates and PRS
#-----------------------------------------

cgtList = c('fit','nedu','tdi') 
for (run_type in cgtList) {
    message('Cognitive Test Type : ',run_type)

    outDir = '/hps/nobackup/dunham/ai-uk-can/analysis/results/prs_interaction'
    outGlm = paste(outDir,paste0(run_type,'_glm_out.txt',sep=""),sep='/')
    outRsq = paste(outDir,paste0(run_type,'_rsq_out.txt',sep=""),sep='/')

    frame.miss.pli.cov = mergeCov(frame.miss.pli,run_type)
    prs.df.slice = getPRScut(prsFile)
    frame.miss.pli.cov.prs.merge = dplyr::left_join(frame.miss.pli.cov,prs.df.slice,by=join_by(f.eid_samples == FID))
    if(run_type != 'tdi'){
        frame.miss.pli.cov.prs = frame.miss.pli.cov.prs.merge[,-c(1,31,32,33,35)] #Drop meta information
    }else{
        frame.miss.pli.cov.prs = frame.miss.pli.cov.prs.merge[,-c(1,35,36,37,39)] #Drop meta information
    
    }

    frame.miss.pli.cov.prs$PRS[is.na(frame.miss.pli.cov.prs$PRS)] = 0
    frame.miss.pli.cov.prs$PRS = scale(frame.miss.pli.cov.prs$PRS)

    #--------------------
    # PTVs, Missense and PRS
    #--------------------

    #PTV, Missense and PRS
    ptv.miss.prs = frame.miss.pli.cov.prs
    ptv.miss.prs.fm = getRegression(ptv.miss.prs,run_type,'glm','ALL')
    
    res.df = data.frame(summary(ptv.miss.prs.fm)$coefficients)
    res.df.2 = data.frame(cbind(rownames(res.df),res.df))
    colnames(res.df.2) = c("Predictors","GLM-Est","StdError","TValue","PValue")
    ptv.miss.prs.res = subset(res.df.2[order(res.df.2$PValue,decreasing=F),],grepl('frame|miss|PRS',Predictors))

    #Only PTV and PRS
    ptv.prs = frame.miss.pli.cov.prs[,!grepl('miss',colnames(frame.miss.pli.cov.prs))]
    ptv.prs.fm = getRegression(ptv.prs,run_type,'glm','ALL')

    #Only Missense and PRS
    miss.prs = frame.miss.pli.cov.prs[,!grepl('frame',colnames(frame.miss.pli.cov.prs))]
    miss.prs.fm = getRegression(miss.prs,run_type,'glm','ALL')

    #Only PRS
    prs = frame.miss.pli.cov.prs[,!grepl('frame|miss',colnames(frame.miss.pli.cov.prs))]
    prs.fm = getRegression(prs,run_type,'glm','ALL')

    #Compute partial R-squared
    reg.gene.df = ptv.miss.prs
    ptv.prs.rsq = rsq.partial(ptv.miss.prs.fm,miss.prs.fm,adj=TRUE,type=c('sse'))$partial.rsq
    miss.prs.rsq = rsq.partial(ptv.miss.prs.fm,ptv.prs.fm,adj=TRUE,type=c('sse'))$partial.rsq
    ptv.miss.rsq = rsq.partial(ptv.miss.prs.fm,prs.fm,adj=TRUE,type=c('sse'))$partial.rsq

    out.prs.rsq = data.frame(c('ptv.prs.rsq','miss.prs.rsq','ptv.miss.rsq'),c(ptv.prs.rsq,miss.prs.rsq,ptv.miss.rsq))
    colnames(out.prs.rsq) = c('Partial_RSquare_Comb','Score')
    
    #------------------------------------------
    # PRS interaction with PTVs and Missense
    #------------------------------------------

    ptv.miss.prs.int = ptv.miss.prs
    ptv.miss.prs.int$ptv_prs = ptv.miss.prs$gene.sum.frame * ptv.miss.prs$PRS
    ptv.miss.prs.int$miss_prs = ptv.miss.prs$gene.sum.miss * ptv.miss.prs$PRS

    ptv.miss.prs.int.fm = getRegression(ptv.miss.prs.int,run_type,'glm','ALL')

    res.df = data.frame(summary(ptv.miss.prs.int.fm)$coefficients)
    res.df.2 = data.frame(cbind(rownames(res.df),res.df))
    colnames(res.df.2) = c("Predictors","GLM-Est","StdError","TValue","PValue")
    ptv.miss.prs.int.res = subset(res.df.2[order(res.df.2$PValue,decreasing=F),],grepl('frame|miss|PRS|prs',Predictors))

    # Printing output table
    out.glm.comb = rbind(ptv.miss.prs.res,ptv.miss.prs.int.res)
    write.table(out.glm.comb,sep='\t',file=outGlm,row.names=F,col.names=T,quote=F)
    write.table(out.prs.rsq,sep='\t',file=outRsq,row.names=F,col.names=T,quote=F)

}

stop()
#------------------------------------------
# PRS interaction with PTV genes
#-----------------------------------------

load('/hps/nobackup/dunham/ai-uk-can/analysis/qualf/mat/cgt.qualf.mat.frame_impact')

gene.ptv.prs.merge = dplyr::left_join(cgt.rm.df,prs.df.slice,by=join_by(f.eid_samples == FID))
gene.ptv.prs = gene.ptv.prs.merge[,-c(1,262:264,266)]
gene.ptv.prs$PRS = scale(gene.ptv.prs$PRS)

gene.ptv.meta = gene.ptv.prs[,c(234:261)]
gene.ptv.meta.fm = getRegression(gene.ptv.meta,'nedu','glm','ALL')
ptv.res = residuals(gene.ptv.meta.fm)
gene.ptv.meta$nedu_resid[is.na(gene.ptv.meta$nedu)==FALSE] = ptv.res

#gene.ptv.meta.res = cbind(gene.ptv.prs[,c(1:233)],gene.ptv.meta$nedu_resid)
gene.ptv.meta.res = gene.ptv.prs #cbind(gene.ptv.prs[,c(1:233)],gene.ptv.meta$nedu_resid)
colnames(gene.ptv.meta.res)[234] = 'nedu'

resList = glmnetcrFunc(gene.ptv.prs,"normal",'nedu')
res.lnz.glm.sort = glmnet_result(resList)


#---
# Full gene, PRS and Interaction term
#---

#gene.ptv.prs.full = gene.ptv.prs$PRS*gene.ptv.prs[,-c(260:261)]
gene.ptv.prs.full = gene.ptv.prs$PRS*gene.ptv.prs[,grepl('NF1',colnames(gene.ptv.prs))]
#colnames(gene.ptv.prs.full) = paste0(colnames(gene.ptv.prs.full),'_prs',sep="")
colnames(gene.ptv.prs.full)[1] = 'NF1_prs'

gene.ptv.prs.int = cbind(gene.ptv.prs,gene.ptv.prs.full)

resList = glmnetcrFunc(gene.ptv.prs.int,'normal','nedu')

#--------------------------------------
# PRS interaction with Missense genes
#--------------------------------------

load('/hps/nobackup/dunham/ai-uk-can/analysis/qualf/mat/cgt.qualf.mat.miss_impact')

gene.miss.prs.merge = dplyr::left_join(cgt.rm.df,prs.df.slice,by=join_by(f.eid_samples == FID))
gene.miss.prs = gene.miss.prs.merge[,-c(1,411:413,415)]
gene.miss.prs$PRS = scale(gene.miss.prs$PRS)

gene.miss.prs.full = gene.miss.prs$PRS*gene.miss.prs[,grepl('ACTL6A|DPP6',colnames(gene.miss.prs))]
colnames(gene.miss.prs.full) = c('ACTL6A_prs','DPP6_prs')

gene.miss.prs.int = cbind(gene.miss.prs,gene.miss.prs.full)



stop()

caseList = c('case1','case2','case3','case4')
#caseList = c('case2')#,'case2','case3','case4')

glmList = list()
glm.df = data.frame(matrix(ncol=10,nrow=0))
colnames(glm.df) = c('run_type','case_id','levels','beta','se','pval','ptv_25','ptv_97','miss_25','miss_97')

count = 1
for (run_type in c('fit','nedu','tdi')){
#for (run_type in c('nedu')){
    
    frame.miss.pli.cov = mergeCov(frame.miss.pli,run_type)
    if(run_type == 'nedu'){
        frame.miss.pli.cov$nedu = frame.miss.pli.cov$nedu
        #frame.miss.pli.cov = frame.miss.pli.cov[,-6]
    }

    prs.df.slice = getPRScut(prsFile)
    levelList = levels(prs.df.slice$SLICE)

    for (case_id in caseList) {
        
        message('-- Processing case IDs: ',case_id)
        
        if(case_id == 'case1'){
            frame.ca = subset(frame.miss.pli.cov,gene.sum.frame==1 & gene.sum.miss ==1)
        }else if (case_id == 'case2'){
            frame.ca = subset(frame.miss.pli.cov,gene.sum.frame==1 & gene.sum.miss ==0)
        }else if (case_id == 'case3'){
            frame.ca = subset(frame.miss.pli.cov,gene.sum.frame==0 & gene.sum.miss ==1)
        }else if (case_id == 'case4'){
            frame.ca = subset(frame.miss.pli.cov,gene.sum.frame==0 & gene.sum.miss ==0)
        }
        
        for (levels in levelList){

            prs.qntl = subset(prs.df.slice,SLICE==levels)
            prs.ref = subset(prs.df.slice,SLICE==3)
            frame.miss.pli.nc = subset(frame.miss.pli.cov,gene.sum.frame ==0 & gene.sum.miss == 0)

            frame.qntl.ca = dplyr::inner_join(prs.qntl,frame.ca,by = join_by(FID == f.eid_samples))
            prs.ref.nc = dplyr::inner_join(prs.ref,frame.miss.pli.nc,by = join_by(FID == f.eid_samples))
           
            if (case_id == 'case4'){
                frame.qntl.ca$gene.sum.frame[frame.qntl.ca$gene.sum.frame==0] = 1
                frame.qntl.ca$gene.sum.miss[frame.qntl.ca$gene.sum.sum==0] = 1
            }
            
            frame.ca.ref.nc = rbind(frame.qntl.ca,prs.ref.nc)
            fm.cgt.cov = frame.ca.ref.nc 
        
            if (run_type == 'tdi') {
                #fm.cgt.cov$tdi = scale(fm.cgt.cov$tdi)
                cgt.rm.df.1 = fm.cgt.cov[,-c(1:5,6)]
                
            }else if (run_type == 'fit'){
                fm.cgt.cov$fit = scale(fm.cgt.cov$fit)
                cgt.rm.df.1 = fm.cgt.cov[,-c(1:5,6)]
                #cgt.rm.df.1 = fm.cgt.cov[,c(7,8,9)]
            }else if (run_type == 'nedu'){
                fm.cgt.cov$nedu = scale(fm.cgt.cov$nedu)
                cgt.rm.df.1 = fm.cgt.cov[,-c(1:5,6)]
            }
            
            reg.rm.df = cgt.rm.df.1
            cat('Processing Levels: ',levels,'\t',dim(reg.rm.df)[1],'\t',dim(reg.rm.df)[2],'\n')
            
            fm = getRegression(reg.rm.df,run_type,"glm","ALL")
            res.df = data.frame(summary(fm)$coefficients)
            res.df.2 = data.frame(cbind(rownames(res.df),res.df))
            colnames(res.df.2) = c("Predictors","GLM-Est","StdError","TValue","PValue")
            res.lnz.glm.sort = res.df.2[order(res.df.2$PValue,decreasing=F),]

            glmList[[run_type]][[case_id]][[levels]] = res.lnz.glm.sort
            beta = subset(res.lnz.glm.sort,grepl('gene.sum',Predictors))[,2]
            se = subset(res.lnz.glm.sort,grepl('gene.sum',Predictors))[,3]
            pval = subset(res.lnz.glm.sort,grepl('gene.sum',Predictors))[,5]
            ptv_conf_int = confint(fm,'gene.sum.frame')
            miss_conf_int = confint(fm,'gene.sum.miss')
            
            glm.df[count,] = c(run_type,case_id,levels,beta,se,pval,ptv_conf_int,miss_conf_int) 
        
            count = count+1
        }
    }

}

write.table(glm.df,sep='\t',file='/hps/nobackup/dunham/ai-uk-can/analysis/results/prs/edu_prs_rv_cgt.txt',row.names=F,col.names=T,quote=F)
