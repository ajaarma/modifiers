###########################################################
#
# Description: Subroutine to compute correlation matrix 
#              for large scale dataset
#
# Author: Ajay A. Kumar (aak@ebi.ac.uk)
#
#
###########################################################

bigcor <- function(x, nblocks = 10, verbose = TRUE, ...){
	
	library(ff, quietly = TRUE)
	NCOL <- ncol(x)

    ## test if ncol(x) %% nblocks gives remainder 0
	if (NCOL %% nblocks != 0) stop("Choose different 'nblocks' so that ncol(x) %% nblocks = 0!")
    
    ## preallocate square matrix of dimension
    ## ncol(x) in 'ff' single format
	corMAT <- ff(vmode = "single", dim = c(NCOL, NCOL))
    
    ## split column numbers into 'nblocks' groups
	SPLIT <- split(1:NCOL, rep(1:nblocks, each = NCOL/nblocks))
    ## create all unique combinations of blocks
	COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
	COMBS <- t(apply(COMBS, 1, sort))
	COMBS <- unique(COMBS)
    
    ## iterate through each block combination, calculate correlation matrix
    ## between blocks and store them in the preallocated matrix on both
    ## symmetric sides of the diagonal
	for (i in 1:nrow(COMBS)) {
        COMB <- COMBS[i, ]
        G1 <- SPLIT[[COMB[1]]]
        G2 <- SPLIT[[COMB[2]]]
        if (verbose) cat("Block", COMB[1], "with Block", COMB[2], "\n")
        flush.console()
        COR <- cor(corMAT[, G1], corMAT[, G2], ...)
        corMAT[G1, G2] <- COR
        corMAT[G2, G1] <- t(COR)
        COR <- NULL
	}
	gc()
	return(corMAT)
}

rmZeroCols <- function(inpDF,run_type){
    
    if (run_type == 'tdi'){
        message('Spliting varDF for: ',run_type)
        rest_ind = c(1:5)
        inpRest = inpDF[,c(1:5)]
        inpVars = inpDF[,-c(1:5)]

    }else{
        rest_ind = c(1)
        inpRest = inpDF[,1]
        inpVars = inpDF[,-1]
    }

    zero.cols = apply(inpVars,2,function(x){
                        #x.na = x[!is.na(x)]
                        tmp = sum(x,na.rm=T)
                        if(tmp >0){
                            return(TRUE)
                        }else{
                            return(FALSE)
                        }
                })

    tmp.df = inpVars[,zero.cols]
    tmp.merge = cbind(inpRest,tmp.df)
    colnames(tmp.merge)[rest_ind] = colnames(inpDF)[rest_ind]
    return(tmp.merge)
}

getCorMat <- function(inpDF){

    inpRest = inpDF[,1]
    inpVars = inpDF[,-1]

    tmpCorMat = mat.or.vec(dim(inpVars)[2],dim(inpVars)[2])
    rownames(tmpCorMat) = colnames(inpVars)
    colnames(tmpCorMat) = colnames(inpVars)
    diag(tmpCorMat) = 1


    tmpJacMat = mat.or.vec(dim(inpVars)[2],dim(inpVars)[2])
    rownames(tmpJacMat) = colnames(inpVars)
    colnames(tmpJacMat) = colnames(inpVars)
    diag(tmpJacMat) = 1

    results = list()

    ncols = dim(tmpCorMat)[1]

    for(i in 1:(ncols-1)){
        print(i)
        k = i+1
        for (j in k:(ncols-1)) {
            r_name = rownames(tmpCorMat)[i]
            c_name = colnames(tmpCorMat)[j]
            cor_val = cor(inpVars[,i],inpVars[,j])
            jac_val = jaccard(inpVars[,i],inpVars[,j])
            
            tmpCorMat[i,j] = round(cor_val,4)
            tmpJacMat[i,j] = round(jac_val,4)
       
        }

    } 
   
    tmpCorMat.full = tmpCorMat+t(tmpCorMat)
    diag(tmpCorMat.full) = diag(tmpCorMat.full)/2
    tmpJacMat.full = tmpJacMat+t(tmpJacMat)
    diag(tmpJacMat.full) = diag(tmpJacMat.full)/2

    results[[1]] = tmpCorMat.full
    results[[2]] = tmpJacMat.full

    return(results)
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


myFunc <- function(query_list){
    
    library("Matrix")
    row.ind = query_list[[1]]
    mat     = query_list[[2]]
    X.df.full = query_list[[3]]
    cluster_num = query_list[[4]]
    out_path = query_list[[5]]

    mat.sub = mat[row.ind,]

    if(class(mat.sub)[1] =='numeric'){
        mat.sub = t(as.matrix(mat.sub))
    }else{
        mat.sub = as.matrix(mat.sub)
    }
    
    print(dim(mat.sub))
    #nrow = dim(query_list[[2]])[1]
   
    consoleFile = paste0(out_path,'/console.txt',collapse="")
    for(i in 1:dim(mat.sub)[1]){
        #cat(cluster_num,'\t',i,"\n",file="/nfs/research/dunham/samples/ukbb/analysis/gmf/cor_mat/console.txt",append=TRUE)
        cat(cluster_num,'\t',i,"\n",file=consoleFile,append=TRUE)

        if(dim(mat.sub)[1]==1){
            row_name = colnames(mat.sub)[row.ind]
        }else{
            row_name = rownames(mat.sub)[i]
        }
        
        col.ind = which(colnames(X.df.full)==row_name)
        row_xdf_vec = as.vector(X.df.full[,col.ind])
        cor_val_vec = apply(X.df.full,2,function(x){
                                cor_val = cor(row_xdf_vec,x)
                                return(cor_val)
                            })
        mat.sub[i,] = cor_val_vec
    }
    return(mat.sub)
}

matrix_init <- function(inpDF, num_node, mode,out_res_dir){#, rowInd=rowInd){

    library("Matrix")
    X.df.rest = inpDF[,1]
    X.df.cov = as.matrix(inpDF[,-1])
    #X.df.cov = as.matrix(inpDF)
    
    mat = X.df.cov

    if(mode==1){
        #mat = sparseMatrix(length(colnames(X.df.cov)),length(colnames(X.df.cov)))
        mat = mat.or.vec(dim(X.df.cov)[2],dim(X.df.cov)[2])
        rownames(mat) = colnames(X.df.cov)
        colnames(mat) = colnames(X.df.cov)
        #mat = mat*1
        nrow = dim(mat)[1]
        ncol = dim(mat)[2]
    }
    
    query.ind = list()
    seq.num = seq(1,nrow,round(nrow/num_node))
    
    for(i in 1:length(seq.num)){
        if(i != length(seq.num)){
            query.ind[[i]] = list(c(seq.num[i]:(seq.num[i+1]-1)),mat,X.df.cov,i,out_res_dir)
        }else {
            query.ind[[i]] = list(c(seq.num[i]:nrow),mat,X.df.cov,i,out_res_dir)
        }
    }
    return(query.ind)
}

make_cluster <- function(query.ind,mode,path){

    save(query.ind,file=paste(path,"/query.ind",sep=""))
    require('snow')
    cl = makeCluster(length(query.ind),type="SOCK")
    #cl = makeCluster(15,type="SOCK")
    cat("after clusater\n")
    tryCatch({
                #print("before cl")
        if(mode==1){
            result = clusterApply(cl,query.ind,myFunc)
        }else{
        print("Inside else, make cluster")
         result = clusterApply(cl,query.ind,myInfoGain)
        }
    },interrupt = function(ex){print(ex)},
        error = function(ex){print(ex)})  
        stopCluster(cl)
        mat.final = result[[1]]
        save(result,file=paste(path,"/CL_result_mat",sep=""))
        for(i in 2:length(query.ind)){
            mat.final = rbind(mat.final,result[[i]])

        }
    return(mat.final)
}

getCovCorMat <- function(inpDF,out_res_dir){

    print(dim(inpDF))
    query.ind = matrix_init(inpDF,20,1,out_res_dir)
    mat.merge = make_cluster(query.ind,1,out_res_dir)
    
    return(mat.merge)
}


