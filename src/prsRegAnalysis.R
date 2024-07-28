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

frame.miss.pli = getPtvMissPli()

#stop()
caseList = c('case1','case2','case3','case4')
#caseList = c('case1')#,'case2','case3','case4')

fmList = list()
glmList = list()
glm.df = data.frame(matrix(ncol=10,nrow=0))
colnames(glm.df) = c('run_type','case_id','levels','beta','se','pval','ptv_25','ptv_97','miss_25','miss_97')

count = 1
for (run_type in c('fit','nedu','tdi')){
#for (run_type in c('tdi')){
    
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
                #cgt.rm.df.1 = fm.cgt.cov[,-c(1:5,6)]
                cgt.rm.df.1 = fm.cgt.cov[,-c(1:4,6)]
                
            }else if (run_type == 'fit'){
                fm.cgt.cov$fit = scale(fm.cgt.cov$fit)
                #cgt.rm.df.1 = fm.cgt.cov[,-c(1:5,6)]
                cgt.rm.df.1 = fm.cgt.cov[,-c(1:4,6)]
                #cgt.rm.df.1 = fm.cgt.cov[,c(7,8,9)]

            }else if (run_type == 'nedu'){
                fm.cgt.cov$nedu = scale(fm.cgt.cov$nedu)
                #cgt.rm.df.1 = fm.cgt.cov[,-c(1:5,6)]
                cgt.rm.df.1 = fm.cgt.cov[,-c(1:4,6)]
            }
            
            reg.rm.df = cgt.rm.df.1
            reg.rm.df$ptv_prs = reg.rm.df$PRS * reg.rm.df$gene.sum.frame
            reg.rm.df$miss_prs = reg.rm.df$PRS * reg.rm.df$gene.sum.miss

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
            
            fmList[[run_type]][[case_id]][[levels]] = fm
            count = count+1
        }
    }

}

#write.table(glm.df,sep='\t',file='/hps/nobackup/dunham/ai-uk-can/analysis/results/prs/edu_prs_rv_cgt.txt',row.names=F,col.names=T,quote=F)
#write.table(glm.df,sep='\t',file='/hps/nobackup/dunham/ai-uk-can/analysis/results/prs/edu_prs_rv_int_cgt.txt',row.names=F,col.names=T,quote=F)
save(fmList,file="/hps/nobackup/dunham/ai-uk-can/analysis/results/prs/edu_prs_rv_int_cgt_fm.list")
