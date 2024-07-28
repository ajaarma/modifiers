###########################################################
#
# Description: Functions to process cognitive Tests in R
#
# Author: Ajay A. Kumar (aak@ebi.ac.uk)
#
#
##########################################################


func.sds.cgt <- function(cgtData){
    "Function to get SDS-CGT vector"

    message('--Processing SDS scores')
    y_ind = c(which(grepl('sds',colnames(cgtData))==TRUE))

    #Digit Substitution Test
    sds_df = data.frame(cgtData)[,y_ind]

    sds.df = data.frame(cbind(cgtData$f.eid_samples,sds_df))
    colnames(sds.df) = c('f.eid_samples','sds')

    Y_sds_vec = data.frame(as.numeric(sds.df[,2]))
    colnames(Y_sds_vec) = c('sds')
    rownames(Y_sds_vec) = as.vector(sds.df$f.eid_samples)
    
    return(Y_sds_vec)

}

func.tmtA.cgt <- function(cgtData){
    "Function to get tmtA-CGT vector"

    message('--Processing tmtA scores')
    y_ind = c(which(grepl('tmta',colnames(cgtData))==TRUE))

    #Qualifications
    tmtA_df = data.frame(cgtData)[,y_ind]

    tmtA.df = data.frame(cbind(cgtData$f.eid_samples,tmtA_df))
    colnames(tmtA.df) = c('f.eid_samples','tmtA')

    Y_tmtA_vec = data.frame(tmtA.df[,2])
    colnames(Y_tmtA_vec) = c('tmtA')
    rownames(Y_tmtA_vec) = as.vector(tmtA.df$f.eid_samples)
    
    return(Y_tmtA_vec)

}

func.tmtB.cgt <- function(cgtData){
    "Function to get tmtA-CGT vector"

    message('--Processing tmtB scores')
    y_ind = c(which(grepl('tmtb',colnames(cgtData))==TRUE))

    #Qualifications
    tmtB_df = data.frame(cgtData)[,y_ind]

    tmtB.df = data.frame(cbind(cgtData$f.eid_samples,tmtB_df))
    colnames(tmtB.df) = c('f.eid_samples','tmtB')

    Y_tmtB_vec = data.frame(tmtB.df[,2])
    colnames(Y_tmtB_vec) = c('tmtB')
    rownames(Y_tmtB_vec) = as.vector(tmtB.df$f.eid_samples)
    
    return(Y_tmtB_vec)

}


func.fit.cgt <- function(cgtData){
    "Function to get FIT-CGT vector"

    message('--Processing FIT scores')

    y_ind = c(which(grepl('fit',colnames(cgtData))==TRUE))

    #Qualifications
    fit_df = data.frame(cgtData)[,y_ind]

    fit.df = data.frame(cbind(cgtData$f.eid,fit_df))
    colnames(fit.df) = c('f.eid_samples','fit')

    Y_fit_vec = data.frame(fit.df[,2])
    colnames(Y_fit_vec) = c('fit')
    rownames(Y_fit_vec) = as.vector(fit.df$f.eid_samples)
    
    return(Y_fit_vec)

}

func.dst.cgt <- function(cgtData){
    "Function to get DST-CGT vector"

    message('--Processing DST score')
    y_ind = c(which(grepl('dst',colnames(cgtData))==TRUE))

    #Qualifications
    dst_df = data.frame(cgtData)[,y_ind]

    dst.df = data.frame(cbind(cgtData$f.eid_samples,dst_df))
    colnames(dst.df) = c('f.eid_samples','dst')

    Y_dst_vec = data.frame(dst.df[,2])
    colnames(Y_dst_vec) = c('dst')
    rownames(Y_dst_vec) = as.vector(dst.df$f.eid_samples)
    
    return(Y_dst_vec)

}
func.rtt.cgt <- function(cgtData){
    "Function to get RTT-CGT vector"

    message('--Processing RTT score')
    y_ind = c(which(grepl('rtt',colnames(cgtData))==TRUE))

    #Qualifications
    rtt_df = data.frame(cgtData)[,y_ind]

    rtt.df = data.frame(cbind(cgtData$f.eid_samples,rtt_df))
    colnames(rtt.df) = c('f.eid_samples',colnames(rtt_df))

    Y_rtt_avg = data.frame(apply(rtt.df,1,function(y){
                                            x = y[-1]
                                            x.tmp = x[!is.na(x)]
                                            if (length(x.tmp)==0){
                                                tmp = NA
                                            }else{
                                                tmp = mean(x.tmp,na.rm=T)
                                            }
                                            return(tmp) 
                                        }
                                 )
                            )

    Y_rtt_vec = data.frame(Y_rtt_avg)
    colnames(Y_rtt_vec) = c('rtt')
    rownames(Y_rtt_vec) = as.vector(rtt.df$f.eid_samples)
    
    return(Y_rtt_vec)

}


func.pmt.cgt <- function(cgtData,isPlot=FALSE){
    "Function to get PMT-CGT vector"

    message('--Processing PMT')

    y_ind = c(which(grepl('pmt',colnames(cgtData))==TRUE))

    #Qualifications
    pmt_df = data.frame(cgtData)[,y_ind]

    pmt.df = data.frame(cbind(cgtData$f.eid_samples,pmt_df))
    colnames(pmt.df) = c('f.eid_samples','pmt_0','pmt_1','pmt_2')

    print(head(pmt.df))
    Y_pmt_sum = data.frame(apply(pmt.df,1,function(y){
                                            x = y[-1]
                                            x.tmp = x[!is.na(x)]
                                            if (length(x.tmp)==0){
                                                tmp = NA
                                            }else{
                                                tmp = sum(x.tmp,na.rm=T)
                                            }
                                            return(tmp) 
                                        }
                                 )
                            )

    Y_pmt_sum = data.frame(Y_pmt_sum)
    colnames(Y_pmt_sum) = c('pmt')
    rownames(Y_pmt_sum) = as.vector(pmt.df$f.eid_samples)
    
    #return(Y_pmt_sum)
    #return(Y_pmt_df)

    if (isPlot){
        Y.pmt.rnds = pmt.df
        Y.pmt.sum.log = log(1+Y_pmt_sum)
        Y.pmt.comb = cbind(Y.pmt.rnds[,-1],Y_pmt_sum,Y.pmt.sum.log)
        colnames(Y.pmt.comb) = c('Round-1','Round-2','Round-3','All-Round-Sum','All-Round-Sum-LogT')
        #hist(Y.pmt.comb,breaks=50)
        return(Y.pmt.comb) 
    }else {
        return(Y_pmt_sum)
    }
}

rmGenesLT10 <- function(cgt.cov.vars.df,gene.cols){

    #cgt.cov.vars.df = cgt.cov.vars.df %>% replace(is.na(.),0)
    a2  = cgt.cov.vars.df[,colnames(cgt.cov.vars.df) %in% gene.cols]
    a2.2 = a2 %>% replace(is.na(.),0)
    a2.sum = apply(a2.2,2,sum)
    a2.sum.10 = a2.sum[a2.sum <10]
    message('--Removing Genes <10 samples :', length(a2.sum.10))
    #message('--Removing Genes <5 samples :', length(a2.sum.10))
    cgt.rm.df.raw = cgt.cov.vars.df[,!colnames(cgt.cov.vars.df) %in% names(a2.sum.10)]

    return(cgt.rm.df.raw)
}

getGeneSums <- function(cgt.rm.df,gene.cols){
    a2 = cgt.rm.df[,colnames(cgt.rm.df) %in% gene.cols]
    a2.2 = a2 %>% replace(is.na(.),0)
    a2.sum = apply(a2.2,2,sum)
    
    a2.sum.df = data.frame(cbind(names(a2.sum),a2.sum))
    colnames(a2.sum.df) = c("geneID","count")
    
    return(a2.sum.df)
}

getGeneSumBurden <- function(cgt.rm.df,gene.cols){
    geneMat = cgt.rm.df[,colnames(cgt.rm.df) %in% gene.cols]
    geneMat.2 = geneMat %>% replace(is.na(.),0)
    geneMat.sum = as.matrix(apply(geneMat.2,1,sum))
    colnames(geneMat.sum) = 'gene.sum'
    
    covMat =  cgt.rm.df[,!colnames(cgt.rm.df) %in% gene.cols]

    cov.geneSum.mat = cbind(covMat,geneMat.sum)
    
    return(cov.geneSum.mat)
}
