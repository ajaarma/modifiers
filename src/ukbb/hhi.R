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


getHhiDF <- function(gmfData,mapData,varMat){

    message('Processing HHI Data')
    #### Get Independent variable House Hold Income (HHI) DF ###

    ukb_var_col = rownames(varMat)
    y_ind = c(which(grepl('^f.738.',colnames(gmfData))==TRUE))

    # Map and sort according to UKB sample ids
    ukb_eid_df = subset(mapData,mapData$V1 %in% ukb_var_col)
    ukb_eid_df2 = ukb_eid_df[match(ukb_var_col,ukb_eid_df$V1),]

    #HHI
    hhi_df = data.frame(gmfData)[,y_ind]

    hhi.df = cbind(gmfData$f.eid,hhi_df)
    colnames(hhi.df)[1] = c('f.eid')

    Y_hhi_df2 = subset(hhi.df,hhi.df$f.eid %in% ukb_eid_df2$V2)
    Y_hhi_df = Y_hhi_df2[match(ukb_eid_df2$V2,Y_hhi_df2$f.eid),]
    Y.hhi.vec = data.frame(Y_hhi_df[,2])

    rownames(Y.hhi.vec) = as.vector(Y_hhi_df[,1])
    colnames(Y.hhi.vec) = "HHI"

    return(Y.hhi.vec)
}

getHHICode <- function(reg.filt.df){

    message('Processing HHI Codes')

    Y.hhi.tmp = reg.filt.df[!is.na(reg.filt.df$HHI),]

    Y.hhi.tmp[Y.hhi.tmp$HHI==1,1] = 'LT-18K'
    Y.hhi.tmp[Y.hhi.tmp$HHI==2,1] = '18K-31K'
    Y.hhi.tmp[Y.hhi.tmp$HHI==3,1] = '31K-52K'
    Y.hhi.tmp[Y.hhi.tmp$HHI==4,1] = '52K-100K'
    Y.hhi.tmp[Y.hhi.tmp$HHI==5,1] = 'GT-100K'

    return(Y.hhi.tmp)
}


