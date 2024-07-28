args=(commandArgs(TRUE))
#set.seed(123)
set.seed(43)
if(Sys.info()['sysname']!='Darwin') {
    loc_4 = "/hps/software/users/dunham/R_lib/4.0.3/"
	require('R.utils',lib=loc_4)
    require('data.table',lib=loc_4)
    require('Matrix',lib=loc_4)
    require('moments',lib=loc_4)
    require('MASS',lib=loc_4)
    require('foreign',lib=loc_4)
    require('withr',lib=loc_4)
    require('lattice',lib=loc_4)
    require('survival',lib=loc_4)
    require('Formula',lib=loc_4)
    require('ggplot2',lib=loc_4)
    require('reshape2',lib=loc_4)
    require('backports',lib=loc_4)
    require('Hmisc',lib=loc_4)
    require('ltm',lib=loc_4)

}else{
	require('R.utils')
	require('data.table')
  	require('Matrix')
  	require('moments')
  	require('MASS')
  	require('foreign')
    require('withr')
    require('lattice')
    require('survival')
    require('Formula')
    require('ggplot2')
    require('backports')
    require('Hmisc')
    require('reshape2')
  require('ltm')
}
sourceDirectory(paste("/homes/aak/scripts/AI-UK/src/ukbb/",sep=""))
####################################################################################
#
# Description: General measures of functioning (GMF) and ID genes
# GMF includes: Qualification( 7 categories); Occupation (9
# categories);Household Income; Townsend Deprivation index
# 
# Independent Variables include: Pathogentic variant carriers, Age, Gender
#
# Regression type: Ordinal regression (MASS)
#
# Author: Ajay Anand Kumar (EMBL-EBI)
#         aak@ebi.ac.uk
#         github.com/ajaarma
####################################################################################


####### General Variable ###############
ukb_path = "/hps/nobackup/dunham/ai-uk-can/analysis"

run_type_list = 'qualf' #args[[1]] #'tmtA' #args[[1]]
test_str = 'syn_impact' #args[[2]] #'frame_impact' #args[[2]]
reg_type = 'glm' #'glmnet'

#UKB GMF data
ukb_cgt_gmf_path = paste0(ukb_path,"/matrices/ukb.cgt.gmf",collapse="")
load(ukb_cgt_gmf_path)

#Load only White-Caucasian samples
ukbDF = subset(ukbDF,f.22006.0.0_ethn==1)

message('Finished loading CGT and GMF data')


#Withdrawn Samples; Remove witjdrawn samples from UKB
wdSamples = read.table('/nfs/research/dunham/samples/ukbb/data/450k/samplesINpheno_ToBeFiltered.txt',header=F)



#Extract col index per CGT-Test type
for (run_type in run_type_list) {
    
    for(case_str in test_str) {
        
        case_id = strsplit(case_str,"_")[[1]][1]
       
        if (grepl('impact',case_str)){
            ukb_var_path = paste0(ukb_path,"/impact_ukbb/",case_id,"/",case_str,"_list",collapse="")
        }else if (grepl('urv',case_str)){
            ukb_var_path = paste0(ukb_path,"/urv_ukbb/",case_id,"/",case_str,"_list",collapse="")
        }

        out_res_dir =  paste0(ukb_path,'/',run_type_list,'/',case_id,collapse='')
         
        sort_sig_path = paste0(out_res_dir,"/",case_str,'_lasso_sort_sig.txt',collapse='')
        glm_obj_path = paste0(out_res_dir,"/",case_str,'_glmnet_lasso.dat',collapse='')
        lassoFitPlot = paste0(out_res_dir,"/",case_str,'_lasso_fit.png',collapse='')
        lambdaMinPlot = paste0(out_res_dir,"/",case_str,'_lambda_min.png',collapse='')
        lamdaOut = paste0(out_res_dir,"/",case_str,'_lambda_values.txt',collapse='')

        ###### Load variants DF ###########
        if (grepl('impact',case_str)){
            ukb_mat_path = paste0(ukb_path,"/impact_ukbb/mat/",case_id,".df",collapse="")
        }else if (grepl('urv',case_str)){
            ukb_mat_path = paste0(ukb_path,"/urv_ukbb/mat/",case_id,".df",collapse="")
        }

        case_mat = load(ukb_mat_path)
        message('Loaded the variant DF: ',run_type,' ',case_str,' ',dim(get(case_mat)))
       
        #Read UKB variants
        load(ukb_var_path)
        case.df = outList[[1]]
        varDF = cbind(rownames(case.df),as.data.frame(as.matrix(case.df)))
        colnames(varDF)[1] = 'f.eid_samples'
        #varDF = rmZeroCols(varDF.cases)
        message('Loaded the variant matrix: ',case_str,' ',dim(varDF))

        message('Filtering out withdrawn samples')
        varDF.filt = varDF[!varDF$f.eid_samples %in% wdSamples$V1,]
        varDF = varDF.filt
        
        gene.cols = colnames(varDF[,-1])
        ###### Remove Withdrawn samples #########

        #### Covariates to be included ###############
        # Gender X-variable

        X_sex_df = getGender(ukbDF) 

        ############ Age as covariate #############
        #X_age_df = getAge(ukbDF)
        X_age_df = getStdAge(ukbDF)

        ########### Age-squared as covariate ########
        X_age_sq = X_age_df^2
        colnames(X_age_sq) = "age_sq"

        ########## Interaction term ##############
        X_sex_age = X_sex_df[,-1] * X_age_df
        colnames(X_sex_age) = "sex_age"
        X_sex_age_sq = X_sex_df[,-1] * X_age_sq
        colnames(X_sex_age_sq) = "sex_age_sq"

        ############ Alcohol Intekae #############
        X_alc_df = getAlcoDF(ukbDF)
        
        #Returns vector with scores: -7,-3,1,2,3,4,5,6,-Inf
        #-Inf: Corresponds to those individuals with all NAs for each
        #instances
         
        ########### Get UKB PCs ###############
        X_pcs_df = getUKBpcs(ukbDF)

        if (run_type == 'tdi'){

            #tmpDF = cbind(cbind(X_sex_df,X_age_df))
            tmpDF = cbind(X_sex_df,X_age_df,X_alc_df,X_age_sq,X_sex_age,X_sex_age_sq,X_pcs_df)
            tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)
            #varMat.cov = dplyr::right_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
            varMat.cov = dplyr::inner_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
            
            message('Extract TDI and covariates related info')
            Y_tdi_vec = getTdiDF(ukbDF,'tdi')
            X_ncars_vec = getTdiDF(ukbDF,'ncars')
            X_nhh_vec = getTdiDF(ukbDF,'nhh')
            X_rent_vec = getTdiDF(ukbDF,'rent')
            X_emp_vec = getTdiDF(ukbDF,'emp')
            
            f.eid_samples = as.character(ukbDF$f.eid_samples)
            tdi.comb.df = as.data.frame(cbind(f.eid_samples,Y_tdi_vec,X_ncars_vec,X_nhh_vec,X_rent_vec,X_emp_vec))
            tdi.varMat.cov.1 = dplyr::inner_join(tdi.comb.df,varMat.cov,by = join_by(f.eid_samples == f.eid_samples))
            
            tdi.varMat.cov = tdi.varMat.cov.1[,c(1:6)]
            tdi.comb.df.2 = subset(tdi.varMat.cov,!is.na(tdi))

            tdi.ncars.outList = func.tdi.ncars.cor(tdi.comb.df.2)
            tdi.ncars.df = tdi.ncars.outList[[1]]

            tdi.nhh.outList.cor = func.tdi.nhh.cor(tdi.comb.df.2)[[2]]
            tdi.nhh.outList = func.tdi.nhh.cor(tdi.ncars.df)
            tdi.nhh.df = tdi.nhh.outList[[1]]
            tdi.nhh.df.cor = tdi.nhh.outList[[2]]; rm(tdi.nhh.outList)

            tdi.rent.outList.cor = func.tdi.rent.cor(tdi.comb.df.2)[[2]]
            tdi.rent.outList = func.tdi.rent.cor(tdi.nhh.df)
            tdi.rent.df = tdi.rent.outList[[1]]; rm(tdi.rent.outList)

            tdi.emp.outList.cor = func.tdi.emp.cor(tdi.comb.df.2)[[2]]
            tdi.emp.outList = func.tdi.emp.cor(tdi.rent.df)
            tdi.emp.df = tdi.emp.outList[[1]]; rm(tdi.emp.outList)

            tdi.varMat.cov.2 = dplyr::left_join(tdi.emp.df,varMat.cov,by = join_by(f.eid_samples == f.eid_samples))

            tdi.cov.vars.df = tdi.varMat.cov.2[,-1]
            rownames(tdi.cov.vars.df) = tdi.varMat.cov.2[,1]
            #rm(tdi.varMat.cov.2); rm(tdi.varMat.cov.1)

            tdi.varMat.cov.2$tdi = scale(tdi.varMat.cov.2$tdi)

            cgt.rm.df = tdi.varMat.cov.2
            save(cgt.rm.df,file=paste0(ukb_path,'/',run_type_list,paste0('/mat/cgt.tdi.mat.',test_str,sep=""),sep=""))
          
            ###### Remmove Genes present in <10 samples ############
            cgt.rm.df.raw = rmGenesLT10(tdi.cov.vars.df,gene.cols)

            cgt.burden.mat = getGeneSumBurden(cgt.rm.df,gene.cols)
            reg.rm.df = cgt.burden.mat[,-1]

            #reg.rm.df = tdi.cov.vars.df
            #reg.rm.df = cgt.rm.df.raw
            reg.rm.df$tdi = scale(reg.rm.df$tdi)

            gene.sum.df = getGeneSums(cgt.rm.df.raw,gene.cols)
            cat('Dimension of reg mat before Reg: ',run_type,' ',case_id,' ',dim(reg.rm.df)[1],' ',dim(reg.rm.df)[2],'\n') 
          
            if (reg_type == 'glm') {
                fm = getRegression(reg.rm.df,"tdi","glm",'ALL')
                res.df = data.frame(summary(fm)$coefficients)
                res.df.2 = data.frame(cbind(rownames(res.df),res.df))
                colnames(res.df.2) = c("Predictors","GLM-Est","StdError","TValue","PValue")
                res.lnz.glm.sort = res.df.2[order(res.df.2$PValue,decreasing=F),]
                res.lnz.glm.sort.full = dplyr::left_join(res.lnz.glm.sort,gene.sum.df,by = join_by(Predictors == geneID))
                write.table(res.lnz.glm.sort.full,sep="\t",file=sort_sig_path,col.names=T,row.names=F,quote=F)
            
            }else if (reg_type == 'glmnet'){

                resList = glmnetcrFunc(reg.rm.df,"normal",run_type)
                resList[[5]] = reg.rm.df
            
                coef.lnz = resList$fitCV$coef.nz
                
                fm = resList[[3]]
                res.df = data.frame(summary(fm)$coefficients)
                res.df.2 = data.frame(cbind(rownames(res.df),res.df))
                colnames(res.df.2) = c("Predictors","GLM-Est","StdError","TValue","PValue")
                
                res.lnz.glm = cbind(coef.lnz,res.df.2[match(coef.lnz$Predictors,res.df.2$Predictors),])
                res.lnz.glm.2 = res.lnz.glm[,-3]
                colnames(res.lnz.glm.2) = c("Predictors","Lasso-Est","GLM-Est","GLM-StdError","TValue","PValue")
                res.lnz.glm.sort = res.lnz.glm.2[order(res.lnz.glm.2$PValue,decreasing=F),]

                ##### Merge Glmnet df to variant df ##########
                res.lnz.glm.sort.full = dplyr::left_join(res.lnz.glm.sort,gene.sum.df,by = join_by(Predictors == geneID))
                write.table(res.lnz.glm.sort.full,sep="\t",file=sort_sig_path,col.names=T,row.names=F,quote=F)

                resList[[4]] = reg.rm.df
                resList[[5]] = res.lnz.glm.sort.full
                save(resList,file=glm_obj_path)
            }
 
        }else if (run_type %in% c("sds","tmtA","tmtB","fit","dst","pmt","rtt","qualf")) {
        
            #X_alc_df = getAlcoDF(cgtData,mapData,varMat)
            #varMat.cov = data.frame(cbind(X_sex_df,X_age_df,X_alc_df,data.frame(varMat)))
            
            if (run_type =="sds"){
               
                Y.sds.vec = func.sds.cgt(ukbDF)
                #tmpDF = cbind(Y.sds.vec,X_sex_df,X_age_df,X_alc_df)
                #tmpDF = cbind(Y.sds.vec,X_sex_df,X_age_df,X_alc_df,X_age_sq,X_pcs_df)
                tmpDF = cbind(Y.sds.vec,X_sex_df,X_age_df,X_alc_df,X_age_sq,X_sex_age,X_sex_age_sq,X_pcs_df)
                
                tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)
                #sds.cov.vars.df = dplyr::left_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
                sds.cov.vars.df = dplyr::inner_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
                #save(cgt.rm.df.raw,file=paste(ukb_path,'/',run_type,'/mat/cgt.',run_type,'.raw.df',sep=''))
                
                sds.cov.vars.df.2 = subset(sds.cov.vars.df, (sds >0 & sds <37)|is.na(sds))
                sds.cov.vars.df.2$sds = scale(sds.cov.vars.df.2$sds)
                
                ###### Remmove Genes present in <10 samples ############
                cgt.rm.df.raw = rmGenesLT10(sds.cov.vars.df.2,gene.cols)

                #save(cgt.rm.df.raw,file=paste(ukb_path,'/',run_type,'/mat/cgt.',run_type,'.raw.df',sep=''))

                cgt.rm.df = cgt.rm.df.raw
                save(cgt.rm.df,file=paste(ukb_path,'/',run_type,'/mat/cgt.',run_type,'.df',sep=''))
                #fm = getRegression(sds.cov.vars.df,"sds","glm","ALL")

            }else if (run_type == 'tmtA'){
                Y.tmtA.vec = func.tmtA.cgt(ukbDF)
                Y.tmtA.log = log(Y.tmtA.vec)
                #tmpDF = cbind(Y.tmtA.log,X_sex_df,X_age_df,X_alc_df)
                tmpDF = cbind(Y.tmtA.log,X_sex_df,X_age_df,X_alc_df,X_age_sq,X_sex_age,X_sex_age_sq,X_pcs_df)

                tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)
                #tmtA.cov.vars.df = dplyr::left_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
                tmtA.cov.vars.df = dplyr::inner_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))

                #tmtA.cov.vars.nna = subset(tmtA.cov.vars.df,!is.na(tmtA))
                tmtA.cov.vars.df$tmtA = scale(tmtA.cov.vars.df$tmtA)
                
                ###### Remmove Genes present in <10 samples ############
                cgt.rm.df.raw = rmGenesLT10(tmtA.cov.vars.df,gene.cols)

                #cgt.rm.df = tmtA.cov.vars.df
                cgt.rm.df = cgt.rm.df.raw
                save(cgt.rm.df,file=paste(ukb_path,'/',run_type,'/mat/cgt.',run_type,'.df',sep=''))
                #fm = getRegression(tmtA.cov.vars.df,"tmtA","glm","ALL")
            
            }else if (run_type == 'tmtB'){
                Y.tmtB.vec = func.tmtB.cgt(ukbDF)
                Y.tmtB.log = log(Y.tmtB.vec)
                #tmpDF = cbind(Y.tmtB.log,X_sex_df,X_age_df,X_alc_df)
                tmpDF = cbind(Y.tmtB.log,X_sex_df,X_age_df,X_alc_df,X_age_sq,X_sex_age,X_sex_age_sq,X_pcs_df)
                tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)
                #tmtB.cov.vars.df = dplyr::left_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
                tmtB.cov.vars.df = dplyr::inner_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
                #tmtB.cov.vars.nna = subset(tmtB.cov.vars.df,!is.na(tmtB))
                tmtB.cov.vars.df$tmtB = scale(tmtB.cov.vars.df$tmtB)
                
                ###### Remmove Genes present in <10 samples ############
                cgt.rm.df.raw = rmGenesLT10(tmtB.cov.vars.df,gene.cols)

                #cgt.rm.df = tmtB.cov.vars.df
                cgt.rm.df = cgt.rm.df.raw
                #fm = getRegression(tmtB.cov.vars.df,"tmtB","glm","ALL")
                save(cgt.rm.df,file=paste(ukb_path,'/',run_type,'/mat/cgt.',run_type,'.df',sep=''))

            }else if (run_type == 'fit'){
                Y.fit.vec = func.fit.cgt(ukbDF)
                #tmpDF = cbind(Y.fit.vec,X_sex_df,X_age_df,X_alc_df)
                tmpDF = cbind(Y.fit.vec,X_sex_df,X_age_df,X_alc_df,X_age_sq,X_sex_age,X_sex_age_sq,X_pcs_df)
                tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)

                fit.cov.df = tmpDF
                fit.cov.df$fit = scale(fit.cov.df$fit)
                save(fit.cov.df, file=paste0('/hps/nobackup/dunham/ai-uk-can/analysis/fit/mat/fit.cov.mat',sep=""))
                
                #fit.cov.vars.df.1 = dplyr::left_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
                fit.cov.vars.df.1 = dplyr::inner_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
                fit.cov.vars.df = fit.cov.vars.df.1
                fit.cov.vars.df$fit = scale(fit.cov.vars.df$fit)

                ###### Remmove Genes present in <10 samples ############
                cgt.rm.df.raw = rmGenesLT10(fit.cov.vars.df,gene.cols)

                cgt.rm.df = cgt.rm.df.raw
                save(cgt.rm.df, file=paste0('/hps/nobackup/dunham/ai-uk-can/analysis/fit/mat/cgt.fit.mat.',test_str,sep=""))
                #fm = getRegression(fit.cov.vars.df,"fit","glm","ALL")
            
            }else if (run_type == 'dst'){
                Y.dst.vec = func.dst.cgt(ukbDF)
                #tmpDF = cbind(Y.dst.vec,X_sex_df,X_age_df,X_alc_df)
                tmpDF = cbind(Y.dst.vec,X_sex_df,X_age_df,X_alc_df,X_age_sq,X_sex_age,X_sex_age_sq,X_pcs_df)
                tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)
                #dst.cov.vars.df.1 = dplyr::left_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
                dst.cov.vars.df.1 = dplyr::inner_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
                dst.cov.vars.df = dst.cov.vars.df.1
                dst.cov.vars.df$dst = scale(dst.cov.vars.df.1$dst)

                ###### Remmove Genes present in <10 samples ############
                cgt.rm.df.raw = rmGenesLT10(dst.cov.vars.df,gene.cols)
                cgt.rm.df = cgt.rm.df.raw
 
                save(cgt.rm.df, file='/hps/nobackup/dunham/ai-uk-can/analysis/dst/mat/cgt.dst.df')
                #fm = getRegression(dst.cov.vars.df,"dst","glm","ALL")

            }else if (run_type == 'pmt'){
                Y.pmt.vec = func.pmt.cgt(ukbDF)
                Y.pmt.log = log(1+Y.pmt.vec)
                #tmpDF = cbind(Y.pmt.log,X_sex_df,X_age_df,X_alc_df)
                tmpDF = cbind(Y.pmt.log,X_sex_df,X_age_df,X_alc_df,X_age_sq,X_sex_age,X_sex_age_sq,X_pcs_df)
                tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)
                #pmt.cov.vars.df = dplyr::left_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
                pmt.cov.vars.df = dplyr::inner_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
                pmt.cov.vars.nna.nz = subset(pmt.cov.vars.df,pmt !=0 | is.na(pmt))

                ###### Remmove Genes present in <10 samples ############
                cgt.rm.df.raw = rmGenesLT10(pmt.cov.vars.nna.nz,gene.cols)
                cgt.rm.df = cgt.rm.df.raw
 
                cgt.rm.df$pmt = scale(cgt.rm.df$pmt)
                save(cgt.rm.df, file='/hps/nobackup/dunham/ai-uk-can/analysis/pmt/mat/cgt.pmt.df')
                #fm = getRegression(pmt.cov.vars.df,"pmt","glm","ALL")
            
            }else if (run_type == 'rtt'){
                Y.rt.vec = func.rtt.cgt(ukbDF)
                #tmpDF = cbind(Y.rt.vec,X_sex_df,X_age_df,X_alc_df)
                tmpDF = cbind(Y.rt.vec,X_sex_df,X_age_df,X_alc_df,X_age_sq,X_sex_age,X_sex_age_sq,X_pcs_df)
                tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)
                #rtt.cov.vars.df.1 = dplyr::left_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
                rtt.cov.vars.df.1 = dplyr::inner_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
                rtt.cov.vars.df.2 = subset(rtt.cov.vars.df.1,rtt>200 & rtt <1200)

                ###### Remmove Genes present in <10 samples ############
                cgt.rm.df.raw = rmGenesLT10(rtt.cov.vars.df.2,gene.cols)
                cgt.rm.df = cgt.rm.df.raw
 
                cgt.rm.df$rtt = scale(log(cgt.rm.df$rtt))
                save(cgt.rm.df,file='/hps/nobackup/dunham/ai-uk-can/analysis/rtt/mat/cgt.rtt.df')

            }else if(run_type == "qualf"){
                
                ukbDF$f.eid_samples = as.character(ukbDF$f.eid_samples)
                Y.qual.vec = getQualDF(ukbDF)
                Y.qual.vec.n1 = subset(Y.qual.vec,  qualf !=-3 )
                Y.qual.edu = mapEduYears(Y.qual.vec.n1)
                #Y.qual.edu = subset(Y.qual.vec.n1,nedu <=22 & nedu >13)
                Y.qual.edu$f.eid_samples = as.character(Y.qual.edu$f.eid_samples)                

                tmpDF = cbind(X_sex_df,X_age_df,X_alc_df,X_age_sq,X_sex_age,X_sex_age_sq,X_pcs_df)
                tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)

                tmpDF.qual = dplyr::inner_join(tmpDF,Y.qual.edu,by=join_by(f.eid_samples == f.eid_samples))
               
                qual.cov.df = tmpDF.qual
                #qual.cov.df$nedu = scale(qual.cov.df$nedu)
                #save(qual.cov.df,file=paste0('/hps/nobackup/dunham/ai-uk-can/analysis/qualf/mat/qualf.cov.mat',sep=""))
                
                edu.cov.vars.df.1 = dplyr::inner_join(varDF,tmpDF.qual,by = join_by(f.eid_samples == f.eid_samples))
                
                edu.cov.vars.df.1$nedu = scale(edu.cov.vars.df.1$nedu)
                edu.cov.vars.df = edu.cov.vars.df.1[,!colnames(edu.cov.vars.df.1) %in% c('qualf')]

                ###### Remmove Genes present in <10 samples ############
                cgt.rm.df.raw = rmGenesLT10(edu.cov.vars.df,gene.cols)
                cgt.rm.df = cgt.rm.df.raw
 
                run_type = 'nedu'

                save(cgt.rm.df,file=paste0('/hps/nobackup/dunham/ai-uk-can/analysis/qualf/mat/cgt.qualf.mat.',test_str,sep=""))

            }
            
            cgt.burden.mat = getGeneSumBurden(cgt.rm.df,gene.cols)
           
            #reg.rm.df = rmZeroCols(cgt.rm.df,run_type)
            #reg.rm.df = cgt.rm.df[,-c(1,202)] #frameshift
            #reg.rm.df = cgt.rm.df[,-1] #frameshift
            reg.rm.df = cgt.burden.mat[,-1]
            gene.sum.df = getGeneSums(cgt.rm.df.raw,gene.cols)
            stop()

            cat('Dimension of reg mat before Reg: ',run_type,' ',case_id,' ',dim(reg.rm.df)[1],' ',dim(reg.rm.df)[2],'\n') 

            if (reg_type == 'glm') {
                
                fm = getRegression(reg.rm.df,run_type,"glm",'ALL')
                res.df = data.frame(summary(fm)$coefficients)
                res.df.2 = data.frame(cbind(rownames(res.df),res.df))
                colnames(res.df.2) = c("Predictors","GLM-Est","StdError","TValue","PValue")
                res.lnz.glm.sort = res.df.2[order(res.df.2$PValue,decreasing=F),]
                res.lnz.glm.sort.full = dplyr::left_join(res.lnz.glm.sort,gene.sum.df,by = join_by(Predictors == geneID))
                write.table(res.lnz.glm.sort.full,sep="\t",file=sort_sig_path,col.names=T,row.names=F,quote=F)
            
            }else if (reg_type == 'glmnet'){

                resList = glmnetcrFunc(reg.rm.df,"normal",run_type)
                resList[[5]] = reg.rm.df
            
                coef.lnz = resList$fitCV$coef.nz
                
                fm = resList[[3]]
                res.df = data.frame(summary(fm)$coefficients)
                res.df.2 = data.frame(cbind(rownames(res.df),res.df))
                colnames(res.df.2) = c("Predictors","GLM-Est","StdError","TValue","PValue")
                
                res.lnz.glm = cbind(coef.lnz,res.df.2[match(coef.lnz$Predictors,res.df.2$Predictors),])
                res.lnz.glm.2 = res.lnz.glm[,-3]
                colnames(res.lnz.glm.2) = c("Predictors","Lasso-Est","GLM-Est","GLM-StdError","TValue","PValue")
                res.lnz.glm.sort = res.lnz.glm.2[order(res.lnz.glm.2$PValue,decreasing=F),]

                ##### Merge Glmnet df to variant df ##########
                res.lnz.glm.sort.full = dplyr::left_join(res.lnz.glm.sort,gene.sum.df,by = join_by(Predictors == geneID))
                write.table(res.lnz.glm.sort.full,sep="\t",file=sort_sig_path,col.names=T,row.names=F,quote=F)

                resList[[4]] = reg.rm.df
                resList[[5]] = res.lnz.glm.sort.full
                save(resList,file=glm_obj_path)
            }
            
        }else if(run_type %in% c("occp","hhi")) {
            
            #X_alc_df = getAlcoDF(cgtData,mapData,varMat)
            #varMat.cov = data.frame(cbind(X_sex_df,X_age_df,X_alc_df,data.frame(varMat)))
            #varMat.cov = data.frame(cbind(X_sex_df,X_age_df,data.frame(varMat)))

            if (run_type == 'qualf'){
                Y.qual.vec = getQualDF(ukbDF)
                Y.qual.edu = mapEduYears(Y.qual.vec)
                stop()
                colnames(Y.qual.vec)[2] = "qualf"

                tmpDF = cbind(Y.qual.vec,X_sex_df,X_age_df,X_alc_df)
                tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)

                tmpDF = tmpDF[,-3]
                qualf.cov.vars.df = dplyr::left_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
                qualf.cov.vars.df.nna = subset(qualf.cov.vars.df, !is.na(qualf))
                reg.filt.df = subset(qualf.cov.vars.df.nna, qualf !=-3 & qualf !=6)

                #Assign Qualification code
                #reg.filt.df.1 = reg.filt.df[,-1]
                reg.fc.df = data.frame(getQualCode(reg.filt.df))
                #reg.fc.df = data.frame(getBinaryQualCode(reg.filt.df))
                #reg.fc.df$Qualifications = factor(reg.fc.df$Qualifications)
                reg.fc.df$qualf = ordered(reg.fc.df$qualf,
                                                   levels=c("A-Level","College","CSE",
                                                   #levels=c("College","School",
                                                            "NVQ","O-Level","Other","NOTA"
                                                            #"NVQ","O-Level","NOTA"
                                                            )
                                                   )
                reg.rm.df = reg.fc.df[,-1]
                stop()
                #resList = glmnetcrFunc(reg.rm.df,"ordinal",run_type)

            }else if(run_type =="occp"){
                Y.occup.vec = getOccupDF(gmfData,mapData,varMat)
                reg.raw.df.tmp = cbind(Y.occup.vec,varMat.cov)
                reg.raw.df = reg.raw.df.tmp[,-which(colnames(reg.raw.df.tmp)=='f.eid')]
                colnames(reg.raw.df)[1] = 'Occupation'

                #Filter the regression DF for Occupation code and -Inf
                reg.filt.1 = reg.raw.df[!is.infinite(reg.raw.df$Occupation),]
                reg.filt.2 = data.frame(getOccupCode(reg.filt.1))
                reg.fc.df = reg.filt.2

                #reg.fc.df$Occupation = factor(reg.fc.df$Occupation)
                reg.fc.df$Occupation = ordered(reg.fc.df$Occupation,
                                                   levels=c("SRManager","ProfOccup",
                                                            "AssocProfTechOccup","AdminOccup",
                                                            "SkillTradeOccup","SelfServiceOccup",
                                                            "SalesCustomOccup","ProcPlantMachineOpt",
                                                            "ElementaryOccup"
                                                            )
                                                   )
            }else if (run_type =='hhi'){
                Y.hhi.vec = getHhiDF(gmfData,mapData,varMat)
                reg.raw.df.tmp = cbind(Y.hhi.vec,varMat.cov)
                reg.raw.df = reg.raw.df.tmp[,-which(colnames(reg.raw.df.tmp)=='f.eid')]
                colnames(reg.raw.df)[1] = 'HHI'

                #Filter the regression DF for HHI and -Inf
                reg.filt.1 = reg.raw.df[reg.raw.df$HHI !=-1,]
                reg.filt.2 = reg.filt.1[reg.filt.1$HHI !=-3,]
                reg.filt.3 = data.frame(getHHICode(reg.filt.2))
                reg.fc.df = reg.filt.3

                reg.fc.df$HHI = factor(reg.fc.df$HHI)
                reg.fc.df$HHI = ordered(reg.fc.df$HHI,
                                           levels=c("LT-18K","18K-31K","31K-52K",
                                                    "52K-100K","GT-100K"
                                                    )
                                        )
            }
         
            reg.rm.df = reg.fc.df[,-1]
            resList = glmnetcrFunc(reg.rm.df,"ordinal",run_type)
            save(resList,file=glm_obj_path)
            #load(glm_obj_path)
           
            s_fit = resList$bic
            s_fit_aic = resList$hat$AIC
            s_fit_bic = resList$hat$BIC
            s_fit_class = table(resList$hat$class)
            df_len = length(resList$nz.coef$beta)
            dfOut = data.frame(names(resList$nz.coef$beta),
                               resList$nz.coef$beta,
                               rep(s_fit,df_len),
                               rep(s_fit_aic,df_len),
                               rep(s_fit_bic,df_len),
                               rep(paste0(names(sort(s_fit_class,decreasing=T)),collapse="_"),df_len)
                               )
            colnames(dfOut) = c('Predictors','OrdinalLassoEstimate','OptLambda','OptAIC','OptBIC','OptClass')
            write.table(dfOut,sep='\t',file=sort_sig_path,col.names=T,row.names=F,quote=F)
        
            
        }
    }
}
