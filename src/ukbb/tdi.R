.libPaths("/hps/software/users/dunham/R_lib/4.0.3/")
require('data.table')
require('Matrix')

########################################################################
#
#Description: Functions to process General Measures of Functioning (GMF)
#           Townsend Deprivation Index (TDI)
#
# Author: Ajay A. Kumar (aak@ebi.ac.uk)
#
#
########################################################################


getTdiDF <- function(gmfData,f_type){

    ######## Get Independent variable House Hold Income (HHI) DF ########

    message('Processing TDI index: ',f_type)

    if (f_type == 'tdi'){
        y_ind = c(which(grepl('^f.189.',colnames(gmfData))==TRUE))
    
    }else if (f_type =='ncars'){

        #X -covariates related to TDI; Household having num of cars
        y_ind = c(which(grepl('^f.728.',colnames(gmfData))==TRUE))
    
    }else if (f_type =='nhh'){
    
        #X -covariates related to TDI; Household having num of house-holds (hh)
        y_ind = c(which(grepl('^f.709.',colnames(gmfData))==TRUE))
    }else if(f_type == 'rent'){

        #X -covariates related to TDI; Own or Rent Accomodation (rent)
        y_ind = c(which(grepl('^f.680.',colnames(gmfData))==TRUE))
    
    }else if(f_type == 'emp'){

        #X -covariates related to TDI; Employment status (emp)
        y_ind = c(which(grepl('^f.6142.',colnames(gmfData))==TRUE))
    
    }
        
    y_cols = colnames(gmfData)[y_ind]

    #TDI
    tmp_df = data.frame(gmfData)[,y_ind]

    tmp.df = data.frame(cbind(gmfData$f.eid_samples,tmp_df))
    colnames(tmp.df) = c('f.eid_samples',y_cols)

    Y_tmp_df = tmp.df
   
    #names(Y.tmp.vec) = as.vector(Y_tmp_df[,1])
    if (f_type =='ncars'){
        Y_proc_vec = data.frame(apply(Y_tmp_df,1,function(x) {
                                        x = x[-1]
                                        tmp = max(x,na.rm=T)
                                        if(is.infinite(tmp)){
                                            tmp = 0
                                        }else {
                                            tmp = tmp
                                        }
                                        return(tmp)
                                    })
                                )
        colnames(Y_proc_vec) = c('ncars')
        rownames(Y_proc_vec) = as.vector(Y_tmp_df[,1])
    }else if (f_type == 'nhh'){
        #Extract max number of individuals
        Y_proc_vec = data.frame(apply(Y_tmp_df,1,function(x){
                                        x = x[-1]
                                        tmp = max(x,na.rm=T)
                                        if(is.infinite(tmp)){
                                            tmp=0
                                        }else {
                                            tmp = tmp
                                        }
                                        return(tmp)
                                     })
                                )
        colnames(Y_proc_vec) = c('nhh')
        rownames(Y_proc_vec) = as.vector(Y_tmp_df[,1])
    }else if (f_type == 'rent'){
        
        #Extract rented or own house individuals
        #Y_proc_vec = Y_tmp_df[,1]
        Y_proc_vec = data.frame(apply(Y_tmp_df,1,function(y){
                                        x = y[-1]
                                        x.tmp = x[!is.na(x)]
                                        if (length(x.tmp)==0){
                                            tmp = 0
                                        
                                        }else if (length(x.tmp)==1){
                                            tmp = x.tmp
                                        
                                        }else {
                                            if(all(x.tmp >0)) {
                                                tmp = min(x.tmp)
                                                #cat(y,'\t',x.tmp,'\t',tmp,'\n')
                                            
                                            }else if (any(x.tmp <0)){
                                                x.tmp.pos = x.tmp[!x.tmp <0]
                                                x.tmp.neg = x.tmp[x.tmp<0]
                                                
                                                if (length(x.tmp.pos) >0){
                                                    x.min.pos = min(x.tmp.pos)
                                                    tmp = x.min.pos
                                                
                                                }else {
                                                    tmp = min(x.tmp.neg)
                                                }
                                                #cat(y,'\t',x.tmp,'\t',x.tmp.pos,'\t',tmp,'\n')
                                            }
                                        }
                                        return(tmp)
                                     })
                                )
        colnames(Y_proc_vec) = c('rent')
        rownames(Y_proc_vec) = as.vector(Y_tmp_df[,1]) 
    
    }else if (f_type=='emp') {
        Y_proc_vec = data.frame(apply(Y_tmp_df,1,function(y){
                                        x = y[-1]
                                        x.tmp = x[!is.na(x)]
                                        if (length(x.tmp)==0){
                                            tmp = 0
                                        }else if (length(x.tmp)==1){
                                            tmp = x.tmp
                                        
                                        }else {
                                            if(all(x.tmp >0)) {
                                                tmp = min(x.tmp)
                                                #cat(y,'\t',x.tmp,'\t',tmp,'\n')
                                            }else if (any(x.tmp <0)){
                                                x.tmp.pos = x.tmp[!x.tmp <0]
                                                x.tmp.neg = x.tmp[x.tmp<0]
                                                
                                                if (length(x.tmp.pos) >0){
                                                    x.min.pos = min(x.tmp.pos)
                                                    tmp = x.min.pos
                                                }else {
                                                    tmp = min(x.tmp.neg)
                                                }
                                            }
                                        }
                                        return(tmp)
                                     })
                                )
        colnames(Y_proc_vec) = c('emp')
        rownames(Y_proc_vec) = as.vector(Y_tmp_df[,1]) 
    }else if (f_type == 'tdi'){
   
        Y_proc_vec = data.frame(Y_tmp_df[,2])
        colnames(Y_proc_vec) = c('tdi')
        rownames(Y_proc_vec) = as.vector(Y_tmp_df[,1])
    }

    return(Y_proc_vec)
}



