###############################################################################
#
#Description: Subroutine to extract gender and stratify into three categories
#
#Author: Ajay A. Kumar (EMBL-EBI) (aak@ebi.ac.uk)
#       
#
###############################################################################

getAge <- function(ageData){
    
    #Subroutine to extract gender and stratify into three categories
    message('Reading the Age data')

    y_ind = c(which(grepl('age',colnames(ageData))==TRUE))
    

    # Map and sort according to UKB sample ids

    #Age
    age_df = data.frame(ageData)[,y_ind]
    age.df = data.frame(cbind(ageData$f.eid_samples, age_df))
    colnames(age.df) = c('f.eid_samples',colnames(age_df))

    X.age.vec = data.frame(apply(age.df,1, function(x){
                            tmp = min(x[-1],na.rm=T)
                            return(tmp)
                            }
                            )
                        )
    rownames(X.age.vec) = as.vector(age.df[,1])
    colnames(X.age.vec) = c('age')

    X.age.tmp = X.age.vec
    X.age.lt50 = data.frame(apply(X.age.tmp,1,function(x){
                                      if (x<=50){
                                        return(1)
                                      }else{
                                          return(0)
                                      }
                                    }
                                 )
                            )

    X.age.50_65 = data.frame(apply(X.age.tmp,1,function(x){
                                      if (x>50 & x<=65 ){
                                        return(1)
                                      }else{
                                          return(0)
                                      }
                                    }
                                 )
                            )

    X.age.gt65 = data.frame(apply(X.age.tmp,1,function(x){
                                      if (x>65){
                                        return(1)
                                      }else{
                                          return(0)
                                      }
                                    }
                                 )
                            )

    X.age.df = cbind(X.age.lt50,X.age.50_65,X.age.gt65)
    colnames(X.age.df) = c("age_lt50","age_50_65","age_gt65")
    rownames(X.age.df) = rownames(X.age.tmp)

    return(X.age.df)
}

getStdAge <- function(ukbDF,anal_type = NULL) {

    #Subroutine to extract gender and stratify into three categories
    message('Reading the Age data')

    y_ind = c(which(grepl('age',colnames(ukbDF))==TRUE))

    # Map and sort according to UKB sample ids Age
    age.df = data.frame(ukbDF)[,y_ind]
    X.age.vec = data.frame(apply(age.df,1, function(x){
                            tmp = min(x,na.rm=T)
                            return(tmp)
                            }
                            )
                        )
    colnames(X.age.vec) = "age"
    rownames(X.age.vec) = as.vector(ukbDF$f.eid_samples)

    return(X.age.vec)

}
