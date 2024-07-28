#######################################################################
# Description: Subroutine to get first 20 PCs of UKBB (Field: 22009)
#
# Date: 07.02.2023
# Author: Ajay A. Kumar (aak@ebi.ac.uk)
# 
######################################################################

getUKBpcs <- function(ukbDF){

    y_ind = c(which(grepl('22009',colnames(ukbDF))==TRUE))

    #First 20 pcs
    ukb.pc.df = data.frame(ukbDF)[,c(y_ind)]

    return(ukb.pc.df)

}
