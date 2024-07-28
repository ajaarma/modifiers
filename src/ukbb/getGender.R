#########################################################
# Description: Get the gender data frame
#
# Author: Ajay A. Kumar (aak@ebi.ac.uk)
#
#
########################################################

getGender <- function(cgmfData){

    #Subroutine to extract covariates for non-ID samples
    message('Processing gender information')

    y_ind = c(which(grepl('sex',colnames(cgmfData))==TRUE))


    #Gender
    sex_df = data.frame(cgmfData)[,y_ind]
    #sex.rfm.df = reformGender(sex_df)
    sex.rfm.df = data.frame(sex_df)

    sex.df = data.frame(cbind(cgmfData$f.eid_samples,sex.rfm.df))
    colnames(sex.df) = c('f.eid_samples','Gender')
    #colnames(sex.df) = c('f.eid')

    return(sex.df)

}

reformGender <- function(gender_df){

    gender_df = as.data.frame(gender_df)

    # Reformulate Gender
    gdMat = mat.or.vec(dim(gender_df)[1],2)
    rownames(gdMat) = rownames(gender_df)
    colnames(gdMat) = c("Female","Male")

    tmp_ind_1 = as.vector(which(gender_df[,1]==0))
    gdMat[tmp_ind_1,1] = 1

    tmp_ind_2 = as.vector(which(gender_df[,1]==1))
    gdMat[tmp_ind_2,2] = 1

    return(as.data.frame(gdMat))

}

