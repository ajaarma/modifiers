.libPaths("/hps/software/users/dunham/R_lib/4.0.3/")
require('data.table')
require('Matrix')

############################################################################
#
#Description: Subroutine to extract Alcohol Intake a day before assessment
#
# Author: Ajay A. Kumar (aak@ebi.ac.uk)
#
# 
#
#
########################################################################


getAlcoDF <- function(alcData){

    #### Get Independent variable Qualification DF ###

    y_ind = c(which(grepl('alco',colnames(alcData))==TRUE))

    #Alochol Intake
    alc_df = data.frame(alcData)[,y_ind]

    alc.df = cbind(alcData$f.eid_samples,alc_df)
    colnames(alc.df)[1] = c('f.eid_samples')

    Y.alc.df = data.frame(apply(alc.df,1,function(y){
                                    x = y[-1]
                                    x.na = x[!is.na(x)]
                                    if (length(x.na)==0){
                                        tmp = 0
                                        return(tmp)
                                    }else{
                                        tmp = sum(as.numeric(x.na))
                                        if(tmp >0){
                                            return(1)
                                        }else{
                                            return(0)
                                        }
                                    }
                                }
                            )
                        )

    rownames(Y.alc.df) = alc.df$f.eid_samples
    colnames(Y.alc.df) = "AlcoholIntake"
    #Y.alc.df[Y.alc.df>0,] = 1

    return(Y.alc.df)
}



