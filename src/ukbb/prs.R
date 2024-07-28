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

