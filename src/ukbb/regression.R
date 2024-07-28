##############################################################################
#
# Description: Subroutines to perform regression GLMs and Ordinal regression
#
# Author: Ajay A. Kumar (aak@ebi.ac.uk)
#
#
##############################################################################

getRegression <- function(inpDF,test_type,reg_type,gene_id){

    require('rsq')

    gene_id.grep = paste("^",gene_id,"_",sep="")
    index = which(grepl(gene_id.grep,colnames(inpDF))==T)

    if(reg_type=="glm"){
        if(test_type=='tdi'){
            message(paste0('Inside subroutine:',reg_type,test_type,sep=" "))
            
            #reg.gene.df = data.frame(cbind(reg.gene.cov,reg.gene.vars))
            reg.gene.df = data.frame(inpDF)
            test_type_index = which(grepl(test_type,colnames(reg.gene.df))==TRUE)

            gfm = glm(formula(paste("tdi",paste(colnames(reg.gene.df)[-test_type_index],collapse="+"),sep="~")),data=reg.gene.df)
            #gfm.sq = rsq.partial(gfm,adj=TRUE,type=c('sse'))

            return(gfm)

        }else if(test_type %in% c('sds','tmtA','tmtB','fit','dst','pmt','rtt','nedu')){
            
            message(paste0('Inside subroutine:',reg_type,test_type,sep=" "))
            
            reg.gene.df = data.frame(inpDF)
            #print(colnames(reg.gene.df))
            test_type_index = which(grepl(test_type,colnames(reg.gene.df))==TRUE)
            gfm =glm(formula(paste(test_type,paste(colnames(reg.gene.df)[-test_type_index],collapse="+"),sep="~")),data=reg.gene.df)

            return(gfm)
        }

    }else if(reg_type == "ordinal"){
        if (test_type %in% c("qualf","occup","hhi")){
            message(paste0('Inside subroutine:',reg_type,test_type,sep=" "))
            

            vars.cols.df = inpDF[,-c(1:5)]
            vars.col.names = colnames(vars.cols.df)
            chr.num = apply(data.frame(vars.col.names),1,function(x) {
                                tmp = strsplit(x,"_")[[1]][2]
                                return(tmp)
                                }
                            )
           
            fm_res = list()
            chr.num = 1
            for (chr in unique(chr.num)) {
                chr.grep = paste0("_",chr,"_",collapse="")
                var.ind = grepl(chr.grep,vars.col.names)
                vars.chr.sub = vars.cols.df[,var.ind]

                reg.gene.df = cbind(inpDF[,1:5],vars.chr.sub)

                #fmt_str_lhs = colnames(reg.gene.df)[1]
                #fmt_str_rhs = paste0(colnames(reg.gene.df)[-1],collapse="+")

                #print(fmt_str_lhs)
                #formula_str = paste0(fmt_str_lhs,"~",fmt_str_rhs,collapse="")

                #fm = polr(formula_str,data=reg.gene.df,Hess=TRUE,method="logistic")
                #fm_res[[chr]] = fm

            }
            #return(fm_res)
            return(reg.gene.df)
        }
    }
}


getOrdinalNet <- function(inpRegDF){
    
    #fit2 = ordinalNet(x=X.df,y=Y.df,family="cumulative",link="logit",parallelTerms=TRUE,nonparallelTerms=FALSE)
    X.df = as.matrix(inpRegDF[,-1])
    Y.df = inpRegDF[,1]

}

glmnetcrFunc <- function(inpDF,reg_type,run_type){

    require('glmnetcr')
    #require('glmnet')
    require('selectiveInference')

    #X.df = inpDF[,-1]
    #Y.df = inpDF[,1]

    if(reg_type == 'ordinal'){
        
        run.ind = which(grepl(run_type,colnames(inpDF))==TRUE)
        xcol.ind = which(!grepl(run_type,colnames(inpDF))==TRUE)

        Y.df = inpDF$qualf
        X.df = data.frame(inpDF[,xcol.ind])
        #colnames(Y.df) = run_type

        #df = X.df %>% replace(is.na(.),0)
        #X.df = df
        #X.mat = makeX(X.df,na.impute=T)
        #X.mat = as.matrix(X.df)

        #if(run_type == 'qualf'){
        #    Y.df = inpDF[,172]
        #    X.df = inpDF[,-172]
        #}

        message('Inside Ordinal Lasso')
        fit.lasso = glmnetcr(x=X.df,y=Y.df,method='forward',trace.it=1,alpha=1)
        #fit.lasso = glmnetcr(x=X.df,y=Y.df,method='backward',trace.it=1,alpha=1)
        #fit.lass = ordinalnet(
        message('After fitting lasso')
        
        fit.bic = select.glmnetcr(fit.lasso)
        message('After computing BIC')
        
        fit.coef = coef(fit.lasso,s=fit.bic)
        message('After computing ceof')
        
        nz.coef = nonzero.glmnetcr(fit.lasso,s=fit.bic)
        message('After computing nz-coef')
        
        hat = fitted(fit.lasso,s=fit.bic)
        message('After computing Hat')
        
        res = list()
        res$fit = fit.lasso
        res$bic = fit.bic
        res$coef = fit.coef
        res$nz.coef = nz.coef
        res$hat = hat
        
        return(res)
    }else if(reg_type =='normal'){

        require('glmnet')
        require('dplyr')

        if(run_type %in% c('sds','tmtA','tmtB','fit','dst','pmt','tdi','rtt','nedu')){
          
            message('Inside normal: ',run_type)
            
            #run.ind = which(grepl(run_type,colnames(X.df))==TRUE)
            run.ind = which(grepl(run_type,colnames(inpDF))==TRUE)
            xcol.ind = which(!grepl(run_type,colnames(inpDF))==TRUE)

            Y.df = data.frame(inpDF[,run.ind])
            X.df = data.frame(inpDF[,xcol.ind])
            colnames(Y.df) = run_type

            message('glmnet for Test type: ',run_type)
            df = X.df %>% replace(is.na(.),0)
            X.df = df
            #X.mat = makeX(X.df,na.impute=T)
            X.mat = as.matrix(X.df)
            y.mk = makeX(data.frame(Y.df))
            m = rowSums(y.mk,na.rm=T)
            Y.na = na.replace(y.mk,m)

            fit.lasso = glmnet(x=X.df,y=Y.na,family="gaussian",trace.it=1,alpha=1)
            fit.cv = cv.glmnet(x=X.mat,y=Y.na,family="gaussian",trace.it=1,alpha=1)
            message('After doing glmnet CV')
            
            #Extract Beta's at Lambda.min
            lambda.min = fit.cv$lambda.min
            coef.lmin = coef(fit.cv,s='lambda.min')
            coef.1se = coef(fit.cv,s='lambda.1se')
            message('After extracting Betas and Lambda-min')

            #Sort the Coef-DF according to Estimate values
            coef.df = cbind(rownames(coef.lmin),round(data.frame(matrix(coef.lmin)),5))
            colnames(coef.df) = c("Predictors","Estimate")
            #coef.nz = coef.df[round(coef.df$Estimate,6)!=0,]
            coef.nz = coef.df[coef.df$Estimate !=0,]
            coef.nz.sort = coef.nz[order(coef.nz$Estimate,decreasing=T),]
            message('After sorting the coefficients according to Estimands')

            #Non-Zero (NZ) Predictors  and Coeffs
            X.nz = X.mat[,colnames(X.mat) %in% coef.nz$Predictors]
            Y.na = Y.na

            #Re-perform regression using glm with LASSO-selected variables
            s.glm.df = cbind(Y.na,X.nz)
            colnames(s.glm.df)[1] = run_type
            s.glm.fit = getRegression(s.glm.df,run_type,'glm','ALL')
            message('After doing GLM regression')

            resList = list()
            resList$fitLasso = fit.lasso
            resList$fitCV = fit.cv
            resList$fitCV$coef.nz = coef.nz
            resList$sglm = s.glm.fit

            return(resList)
        }
    }

}

getLassoPvalue <- function(lasso_fit){

    require('selectiveInference')
    set.seed(43)
    n=50
    p=10
    sigma=1

    x = matrix(rnorm(n*p),n,p)

}


getRegFitPlot <- function(resList,lassoFitPlot,lambdaMinPlot){
    #### subroutine to get regression fit and Lmabda min plots #####
    message('Generating regression LASSO and Lambda-Min plots')

    lassoFit = resList[[1]]
    fitCV = resList[[2]]

    png(lassoFitPlot)
    plot(lassoFit)
    dev.off()

    png(lambdaMinPlot)
    plot(fitCV)
    dev.off()
}
