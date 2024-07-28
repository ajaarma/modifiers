##############################################################################
#
#Description: Subroutine to find the multi-collinearity present in the data
#
#Author: Ajay A. Kumar (aak@ebi.ac.uk)
#
#
##############################################################################

getMultCol <- function(inpDF,chrNum){

    library('reshape2')
    pred.df = inpDF[,-1]
    pred.rest.df = pred.df[,c(1:4)]
    pred.vars.df = pred.df[,-c(1:4)]

    chrNum.grep = paste0("_",chrNum,"_",collapse="")
    chr.vars.df = pred.vars.df[,grepl(chrNum.grep,colnames(pred.vars.df))]
   
    chr.cor.mat = cor(chr.vars.df)

    chr.cor.melt = melt(chr.cor.mat)
    #chrList = c(1:22,X,Y)
    #chr.cor.melt.50_01 = chr.cor.melt[chr.cor.melt$value <1 & chr.cor.melt$value >=0.50,]
    #chr.cor.melt.20_01 = chr.cor.melt[chr.cor.melt$value <1 & chr.cor.melt$value >=0.20,]
    chr.cor.melt.05_01 = chr.cor.melt[chr.cor.melt$value <1 & chr.cor.melt$value >=0.05,]

    return(chr.cor.melt)

}

jac_sim <- function(reg.rm.df,var1,var2){

    library('jaccard')

    sim_jac = jaccard(reg.rm.df[,which(colnames(reg.rm.df)==var1)],
                      reg.rm.df[,which(colnames(reg.rm.df)==var2)]
                      )
    sim_cor = cor(reg.rm.df[,which(colnames(reg.rm.df)==var1)],
                      reg.rm.df[,which(colnames(reg.rm.df)==var2)]
                      )

    sim_res = list()
    sim_res[[1]] = sim_jac
    sim_res[[2]] = sim_cor
    return(sim_res)
}

getSimVal <- function(inpDF,var1,var2){

    tmp = inpDF[,which(colnames(inpDF)==var1) & which(colnames(inpDF)==var2)]

    return(tmp)

}
