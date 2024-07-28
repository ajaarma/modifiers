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


getQualDF <- function(gmfData){


    message('--Processing Qualification Data')
    #### Get Independent variable Qualification DF ###
    
    #y_ind = c(which(grepl('^f.6138.',colnames(gmfData))==TRUE))
    y_ind = c(which(grepl('qualf|nedu',colnames(gmfData))==TRUE))
    

    #Qualifications
    qual_df = data.frame(gmfData)[,y_ind]

    qual.df = cbind(gmfData$f.eid_samples,qual_df)
    colnames(qual.df)[1] = c('f.eid_samples')

    Y_qual_df = qual.df[,-c(1,2)]

    Y.qual.vec = data.frame(apply(Y_qual_df,1,function(x){
                            #x1 = x[-1]
                            x1 = x[1]
                            #x1.nna = x1[!is.na(x1)]
                            #if(length(x1.nna) ==0){
                            #    tmp = NA
                            #}else{
                            #    tmp = min(x1.nna)
                            #}
                            #return(tmp)
                            return(x1)
                        }))

    #Y.qual.df = cbind(qual.df$f.eid_samples,Y.qual.vec[,1])
    Y.qual.df = cbind(qual.df[,c(1,2)],Y.qual.vec[,1])
    #rownames(Y.qual.df) = as.vector(Y_qual_df[,1])
    colnames(Y.qual.df) = c("f.eid_samples","nedu","qualf")

    return(Y.qual.df)
}


getQualCode <- function(tmp1){

    message('--Processing Qualification Codes')
    Y.qual.tmp = tmp1
    Y.qual.tmp[Y.qual.tmp$qualf==1,'qualf'] = 'College'
    Y.qual.tmp[Y.qual.tmp$qualf==2,'qualf'] = 'A-Level'
    Y.qual.tmp[Y.qual.tmp$qualf==3,'qualf'] = 'O-Level'
    Y.qual.tmp[Y.qual.tmp$qualf==4,'qualf'] = 'CSE'
    Y.qual.tmp[Y.qual.tmp$qualf==5,'qualf'] = 'NVQ'
    Y.qual.tmp[Y.qual.tmp$qualf==6,'qualf'] = 'Other'
    Y.qual.tmp[Y.qual.tmp$qualf==-7,'qualf'] = 'NOTA'

    return(Y.qual.tmp)
}
 
getBinaryQualCode <- function(tmp1){

    message('--Processing Qualification Codes')
    Y.qual.tmp = tmp1
    Y.qual.tmp[Y.qual.tmp$qualf==1,'qualf'] = 'College'
    Y.qual.tmp[Y.qual.tmp$qualf==2,'qualf'] = 'School'
    Y.qual.tmp[Y.qual.tmp$qualf==3,'qualf'] = 'School'
    Y.qual.tmp[Y.qual.tmp$qualf==4,'qualf'] = 'School'
    Y.qual.tmp[Y.qual.tmp$qualf==5,'qualf'] = 'School'
    Y.qual.tmp[Y.qual.tmp$qualf==6,'qualf'] = 'School'
    Y.qual.tmp[Y.qual.tmp$qualf==-7,'qualf'] = 'School'

    return(Y.qual.tmp)
}
 
mapEduYears <- function(tmp1){

    message('--Mapping Education Years to Qualification Codes')
    #Y.qual.tmp = subset(tmp1,nedu >4 | is.na(nedu))
    #Y.qual.tmp = subset(tmp1,qualf | is.na(nedu))
    Y.qual.tmp = tmp1

    
    ##### Qualf code 1 ##########
    #ind = which(Y.qual.tmp$qualf ==1)
    #Y.qual.tmp[Y.qual.tmp$qualf==1,'nedu'] = 21
    #Y.qual.tmp[ind,'nedu'] = 0
    a1 = subset(Y.qual.tmp,qualf==1)
    a1$edu_years = 20

    ### Qualf code 2
    a2 = subset(Y.qual.tmp,qualf==2)
    #a2$nedu = mean(a2$nedu,na.rm=T) - a2$nedu
    #a2$nedu = 18 - a2$nedu
    a2$edu_years = 13

    #### Qualf code : 3
    a3 = subset(Y.qual.tmp,qualf==3)
    #a3$nedu = mean(a3$nedu,na.rm=T) - a3$nedu
    #a3$nedu = 16 - a3$nedu
    a3$edu_years = 10

    #### Qualf code : 4
    a4 = subset(Y.qual.tmp,qualf==4)
    #a4$nedu = mean(a4$nedu,na.rm=T) - a4$nedu
    #a4$nedu = 16 - a4$nedu
    a4$edu_years = 10

    #### Qualf code: 5
    a5 = subset(Y.qual.tmp,qualf==5)
    #a5$nedu = mean(a5$nedu,na.rm=T) - a5$nedu
    #a5$nedu = 15 - a5$nedu
    a5$edu_years = a5$nedu - 5

    ##### Qualf code: 6
    a6 = subset(Y.qual.tmp,qualf==6)
    a6$edu_years = 15

    #### Qualf code : -7
    a7 = subset(Y.qual.tmp,qualf==-7)
    a7$edu_years = 7
    #a7$nedu = round(mean(a7$nedu,na.rm=T)) - a7$nedu
    #a7$nedu =  - a7$nedu


    #Y.qual.comb = unique(rbind(a2,a3,a4,a5,a7))
    Y.qual.comb = unique(rbind(a1,a2,a3,a4,a5,a6,a7))

    Y.qual.comb = Y.qual.comb[,!colnames(Y.qual.comb) %in% c('nedu')]
    colnames(Y.qual.comb)[3] = 'nedu'

    return(Y.qual.comb)
}
   
 

