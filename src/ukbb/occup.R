.libPaths("/hps/software/users/dunham/R_lib/4.0.3/")
require('data.table')
require('Matrix')

########################################################################
#
#Description: Functions to process General Measures of Functioning (GMF)
#
# Author: Ajay A. Kumar (aak@ebi.ac.uk)
#
#
########################################################################

getOccupDF <- function(gmfData,mapData,varMat){

    message('Processing Occupation data')

    #### Get Independent variable Occupation DF ###
    ukb_var_col = rownames(varMat)
    y_ind = c(which(grepl('^f.22617.',colnames(gmfData))==TRUE))

    # Map and sort according to UKB sample ids
    ukb_eid_df = subset(mapData,mapData$V1 %in% ukb_var_col)
    ukb_eid_df2 = ukb_eid_df[match(ukb_var_col,ukb_eid_df$V1),]

    #Occupations
    occup_df = data.frame(gmfData)[,y_ind]

    occup.df = cbind(gmfData$f.eid,occup_df)
    colnames(occup.df)[1] = c('f.eid')

    Y_occup_df2 = subset(occup.df,occup.df$f.eid %in% ukb_eid_df2$V2)
    Y_occup_df = Y_occup_df2[match(ukb_eid_df2$V2,Y_occup_df2$f.eid),]

    Y.occup.vec = data.frame(apply(Y_occup_df,1,function(y){
                                x = y[-1]
                                x.na = x[!is.na(x)]
                                tmp = min(x.na)
                                return(tmp)
                                }
                                )
                            )

    rownames(Y.occup.vec) = as.vector(Y_occup_df[,1])
    colnames(Y.occup.vec) = c("Occupation")

    return(Y.occup.vec)
}


getOccupCode <- function(reg.filt.df){

    message('Processing Occupation code')
    Y.occup.tmp = reg.filt.df
    Y.occup.tmp[,1] = apply(data.frame(Y.occup.tmp[,1]),1,function(x){
                                tmp = as.numeric(strsplit(as.character(x),"")[[1]][1])
                                return(tmp)
                        })
    colnames(Y.occup.tmp)[1] = 'Occupation'
    Y.occup.tmp[Y.occup.tmp$Occupation==1,1] = 'SRManager'
    Y.occup.tmp[Y.occup.tmp$Occupation==2,1] = 'ProfOccup'
    Y.occup.tmp[Y.occup.tmp$Occupation==3,1] = 'AssocProfTechOccup'
    Y.occup.tmp[Y.occup.tmp$Occupation==4,1] = 'AdminOccup'
    Y.occup.tmp[Y.occup.tmp$Occupation==5,1] = 'SkillTradeOccup'
    Y.occup.tmp[Y.occup.tmp$Occupation==6,1] = 'SelfServiceOccup'
    Y.occup.tmp[Y.occup.tmp$Occupation==7,1] = 'SalesCustomOccup'
    Y.occup.tmp[Y.occup.tmp$Occupation==8,1] = 'ProcPlantMachineOpt'
    Y.occup.tmp[Y.occup.tmp$Occupation==9,1] = 'ElementaryOccup'

    return(Y.occup.tmp)
}
