args=(commandArgs(TRUE))
set.seed(123)
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

getCgtTest <- function(cgtData,mapData,ukb_var_col,x_ind,y_ind){
    # Subroutine to extract Cognitive Test vector corresponding to
    # given UKB sample set


    ukb_eid_df = subset(mapData,mapData$V1 %in% ukb_var_col)
    ukb_eid_df2 = ukb_eid_df[match(ukb_var_col,ukb_eid_df$V1),]

    # Get the CG Test field w.r.t ordered according to sample id
    cgt_df = data.frame(cgtData)[,c(x_ind,y_ind)]
    
    #Employee DF
    cgt_eid_df = as.data.frame(cgt_df[,1])
    colnames(cgt_eid_df) = c('f.eid')

    # cgt gender DF
    cgt_sex_df = reformGender(cgt_df[,2])

    #Alcohol consumption a day before
    cgt_alc_1 = apply(data.frame(cgt_df[,c(3:6)]),1,function(x){y = sum(x,na.rm=T);return (y)})
    cgt_alc_1[cgt_alc_1>0] = 1
    cgt_alc_df = data.frame(cgt_alc_1)
    colnames(cgt_alc_df) = c("AlcoholIntake")

    #CGT Test Matrix
    cgt_pmt_df = cgt_df[,c(7:dim(cgt_df)[2])]

    cgt_X_pmt_tmp_df = cbind(cgt_eid_df,cgt_sex_df,cgt_alc_df,cgt_pmt_df)
    cgt_X_pmt_df2 = subset(cgt_X_pmt_tmp_df,cgt_X_pmt_tmp_df$f.eid %in% ukb_eid_df2$V2)
    cgt_X_pmt_df = cgt_X_pmt_df2[match(ukb_eid_df2$V2,cgt_X_pmt_df2$f.eid),]


    # Get DI and HI
    #hi_df = data.frame(cgtData)[,c(3:7)]
    #hi_df_2 = reformHImat(hi_df)

    out_list = list(ukb_eid_df2,cgt_X_pmt_df)

    return(out_list)
}

reformHImat <- function(hi_df){
    # Reformulate caregorical Household Income data frame as continuous variable

    hiMat = mat.or.vec(dim(hi_df)[1],5)
    rownames(hiMat) = rownames(hi_df)
    colnames(hiMat) = c("lt_18K","bt_18K_31K","bt_31K_52K","bt_52K_100K","gt_100K")

    for (i in c(1:dim(hiMat)[2])) {
        tmp_ind_1 = as.vector(which(hi_df[,1]==i))
        hiMat[tmp_ind_1,i] = 1

        tmp_ind_2 = as.vector(which(hi_df[,2]==i))
        hiMat[tmp_ind_2,i] = 1

        tmp_ind_3 = as.vector(which(hi_df[,3]==i))
        hiMat[tmp_ind_3,i] = 1

        tmp_ind_4 = as.vector(which(hi_df[,4]==i))
        hiMat[tmp_ind_4,i] = 1

    }

    return(as.data.frame(hiMat))

}

getPMTvector <- function(cg_test_df,map_sort_df,run_type){
    # Subroutine to forumate Reponse vector for pairs-matchin test
    # Columns 3-5 indicate number of rounds: 1,2,3
    # Value 0 indicate participant made no mistake
    # Cgnitve score: Sum all the three round score 
    # Cognitive score: For participant no particpating in any round will have NA
    # (missing value)

    cg_test_df = data.frame(cg_test_df)

    na_df = cg_test_df[apply(cg_test_df,1,function(x){y = all(is.na(x[5:dim(cg_test_df)[2]]));return(y)}),]
    na_df_Y = na_df[,c(1,dim(cg_test_df)[2])]
    colnames(na_df_Y) = c("UKB_ID","CGT_Score")

    rest_df = cg_test_df[!apply(data.frame(cg_test_df[,c(5:dim(cg_test_df)[2])]),1,function(x){y = all(is.na(x));return(y)}),]
    rest_df_sum = data.frame(apply(rest_df,1,function(x){y = sum(x[5:dim(cg_test_df)[2]],na.rm=T);return(y)}))
    rest_df_Y.tmp = cbind(rest_df[,1],rest_df_sum)
    colnames(rest_df_Y.tmp) = c("UKB_ID","CGT_Score")

    cgt_df_Y.tmp = rbind(na_df_Y,rest_df_Y.tmp)

    cgt_df_Y = cgt_df_Y.tmp[match(map_sort_df$V2,cgt_df_Y.tmp$UKB_ID),]

    return(cgt_df_Y)
}


getGLM <- function(varMat,cgt_X_pmt_df,cg_test_Y,gene_id,run_type){
    # Subroutine to perform GLM regression for variants of given genes

    gene_id.grep = paste("_",gene_id,"$",sep="")
    index = which(grepl(gene_id.grep,colnames(varMat))==T)

    varGeneMat = data.frame(cbind(varMat[,index],cgt_X_pmt_df[,c(2:4)]))
    colnames(varGeneMat) = c(colnames(varMat)[index],"Female","Male","AlcoholIntake")
    colnames(varGeneMat) = gsub("^","X",colnames(varGeneMat))
   

    if (run_type =='pmt'){
        pmtVarMat = cbind(log(1+cg_test_Y$CGT_Score),varGeneMat)
    }else if(run_type =='tmtA'){
        pmtVarMat = cbind(log(cg_test_Y$CGT_Score),varGeneMat)
    }else if(run_type == 'tmtB'){
        pmtVarMat = cbind(log(cg_test_Y$CGT_Score),varGeneMat)
    }
    else {
        pmtVarMat = cbind(cg_test_Y$CGT_Score,varGeneMat)
    }
        
    colnames(pmtVarMat)[1] = "CGT_SCORE"

    gfm = glm(formula(paste("CGT_SCORE",paste(colnames(varGeneMat),collapse="+"),sep="~")),data=pmtVarMat)

    return(gfm)
}



####### General Variable ###############
ukb_path = "/nfs/research/dunham/samples/ukbb/"
#ukb_path = "/Users/ajaykumar/sshfs/ebi/samples/ukbb"

batchNum = 'merge_all'

run_type_list = 'qualf' #args[[1]]
test_str =  '_vt_urcgte25_hc6' #args[[2]]
id_int_str = 'id' #args[[3]]
gene_id =  'ALL' #args[[4]]
out_res_dir =  paste0('/hps/nobackup/dunham/ukbb/id_genes_comb/',run_type_list,'/',test_str,collapse='')

if (id_int_str=="int"){
    id_int_type = 'int_genes'
    id_int_merge_type = 'int_genes_ukbb_merge_all'
}else if (id_int_str == 'id'){
    id_int_type = 'id_genes'
    id_int_merge_type = 'id_genes_ukbb_merge_all'
}

if (test_str=='_vt_rc25_hc6'){
    gene_id = 'Case1'
}else if(test_str == '_vt_urcgte25_hc6'){
    gene_id = 'Case2'
}else if(test_str == '_vt_urclte25_hc6'){
    gene_id = 'Case3'
}

#out_res_dir = paste0('/hps/nobackup/dunham/ukbb/',id_int_type,'/',run_type_list,"/",test_str,collapse="") #args[[5]]
#system(paste0('mkdir -p ',out_res_dir,collapse=""))

#UKB GMF data
ukb_gmf_path = paste0(ukb_path,"/data/gmf_test/samples_gmf_test.txt",collapse="")
gmfData = fread(ukb_gmf_path,sep="\t",header=T,stringsAsFactor=F,quote="")
#gmfData = data.frame(gmfData)
gmf_col = colnames(gmfData)
message('Finished reading GMF data')

#UKB CGT data
ukb_cgt_path = paste0(ukb_path,"/data/cg_test/samples_cg_test_v3.txt",collapse="")
cgtData = fread(ukb_cgt_path,sep="\t",header=T,stringsAsFactor=F,quote="")
cgt_col = colnames(cgtData)
message('Finished reading CGT data')

#UKB ICD10 Map path
ukb_map_path = paste0(ukb_path,"/data/icd10/",batchNum,"_map.txt",collapse="")
mapData = fread(ukb_map_path,sep="\t",header=F,stringsAsFactor=F,quote="")
#mapData = data.frame(mapData)

#UKB non-ID samples
ukb_ex_sample = paste0(ukb_path,"/data/icd10/NonNeuroF7_samples_50K.txt",collapse="")
ukbNonID = fread(ukb_ex_sample,sep='\t',header=F,stringsAsFactor=F,quote="")
#ukbNonId = data.frame(ukbNonID)
id_samples = ukbNonID$V2


#Extract col index per CGT-Test type
for (run_type in run_type_list) {
    
    for(test_id in test_str) {
        
        ukb_var_path = paste0(ukb_path,"/analysis/filterSNV/",id_int_merge_type,
                              '/merge_all',test_id,".txt",collapse="")
        #sig_path = paste0(out_res_dir,"/",gene_id,'_sig.txt',collapse='')
        sort_sig_path = paste0(out_res_dir,"/",gene_id,'_sort_sig.txt',collapse='')
        fm_obj_path = paste0(out_res_dir,"/",gene_id,'_fm_obj.dat',collapse='')

        #Read UKB variants
        varData = fread(ukb_var_path,sep="\t",header=F,stringsAsFactor=F,quote="")
       
        #Get variant-sample Matrix (rows=UKBB non-ID, col=variants)
        varMat =  getVarMat(varData,mapData,id_samples)
        col.name.vec = apply(data.frame(colnames(varMat)),1,function(x){
                                        strs = strsplit(x,"_")[[1]]
                                        col_str = paste0(strs[3],"_",
                                                         strs[1],"_",
                                                         strs[2],collapse="")
                                        return(col_str)
                                        }
                            )
        colnames(varMat) = col.name.vec

        message('Finished building variant matrix')

        #### Covariates to be included ###############
        # Gender X-variable

        X_sex_df = getGender(cgtData,mapData,varMat) 

        ############ Age as covariate #############
        X_age_df = getAge(cgtData,mapData,varMat)
        

        #Returns vector with scores: -7,-3,1,2,3,4,5,6,-Inf
        #-Inf: Corresponds to those individuals with all NAs for each
        #instances
        if (run_type == 'tdi'){

            varMat.cov = data.frame(cbind(X_sex_df,X_age_df,data.frame(varMat)))

            message('Extract TDI and covariates related info')
            Y_tdi_vec = getTdiDF(gmfData,mapData,varMat,'tdi')
            X_ncars_vec = getTdiDF(gmfData,mapData,varMat,'ncars')
            X_nhh_vec = getTdiDF(gmfData,mapData,varMat,'nhh')
            X_rent_vec = getTdiDF(gmfData,mapData,varMat,'rent')
            X_emp_vec = getTdiDF(gmfData,mapData,varMat,'emp')

            tdi.comb.df = as.data.frame(cbind(Y_tdi_vec,X_ncars_vec,X_nhh_vec,X_rent_vec,X_emp_vec))
            tdi.comb.df.2 = tdi.comb.df[!is.na(tdi.comb.df$tdi),]

            tdi.ncars.outList = func.tdi.ncars.cor(tdi.comb.df.2)
            tdi.ncars.df = tdi.ncars.outList[[1]]

            tdi.nhh.outList = func.tdi.nhh.cor(tdi.comb.df.2)
            tdi.nhh.df = tdi.nhh.outList[[1]]

            tdi.rent.outList = func.tdi.rent.cor(tdi.comb.df.2)
            tdi.rent.df = tdi.rent.outList[[1]]

            tdi.emp.outList = func.tdi.emp.cor(tdi.comb.df.2)
            tdi.emp.df = tdi.emp.outList[[1]]

            message('Intersection between tdi & covariates')
            tdi.cov.common = Reduce(intersect,list(rownames(Y_tdi_vec), rownames(tdi.ncars.df),
                                   rownames(tdi.nhh.df), rownames(tdi.rent.df),
                                   rownames(tdi.emp.df)
                                    )
                             )
            tdi.vec.df2 = Y_tdi_vec[match(tdi.cov.common,rownames(Y_tdi_vec)),]
            tdi.ncars.df2 = tdi.ncars.df[match(tdi.cov.common,rownames(tdi.ncars.df)),]
            tdi.nhh.df2 = tdi.nhh.df[match(tdi.cov.common,rownames(tdi.nhh.df)),]
            tdi.rent.df2 = tdi.rent.df[match(tdi.cov.common,rownames(tdi.rent.df)),]
            tdi.emp.df2 = tdi.emp.df[match(tdi.cov.common,rownames(tdi.emp.df)),]
            tdi.varMat.cov.df2 = varMat.cov[match(tdi.cov.common,varMat.cov$f.eid),]

            tdi.cov.df = data.frame(cbind(tdi.vec.df2,tdi.ncars.df2$ncars,
                                          tdi.nhh.df2$nhh,tdi.rent.df2$rent,
                                          tdi.emp.df2$emp
                                          )
                                    )
            colnames(tdi.cov.df) = c("tdi","ncars","nhh","rent","emp")
            rownames(tdi.cov.df) = tdi.cov.common
            
            tdi.cov.vars.df = cbind(tdi.cov.df,tdi.varMat.cov.df2[,-1])
            
            #fm = getRegression(tdi.cov.vars.df,"tdi","glm",gene_id)
            fm = getRegression(tdi.cov.vars.df,"tdi","glm","ALL")
            save(fm,file=fm_obj_path)
            
            res.df = data.frame(summary(fm)$coefficients)
            res.df[,4] = as.numeric(res.df[,4])
            res.df.sort = res.df[order(res.df[,4],decreasing=F),]
            res.df.sort.2 = data.frame(cbind(rownames(res.df.sort),res.df.sort))
            colnames(res.df.sort.2) = c("Variable","Estimate","StdError","T-Value","P-Value")
            write.table(res.df.sort.2,sep='\t',file=sort_sig_path,col.names=T,row.names=F,quote=F)

        }else if (run_type %in% c("sds","tmtA","tmtB","fit","dst","pmt")) {
        
            X_alc_df = getAlcoDF(cgtData,mapData,varMat)
            varMat.cov = data.frame(cbind(X_sex_df,X_age_df,X_alc_df,data.frame(varMat)))
            
            if (run_type =="sds"){
                
                Y.sds.vec = func.sds.cgt(cgtData,mapData,varMat)
                sds.cov.vars.df = data.frame(cbind(Y.sds.vec,varMat.cov[,!(colnames(varMat.cov) %in% "f.eid")]))
                fm = getRegression(sds.cov.vars.df,"sds","glm","ALL")

            }else if (run_type == 'tmtA'){
                Y.tmtA.vec = func.tmtA.cgt(cgtData,mapData,varMat)
                Y.tmtA.log = log(Y.tmtA.vec)
                tmtA.cov.vars.df = data.frame(cbind(Y.tmtA.log,varMat.cov[,!(colnames(varMat.cov) %in% "f.eid")]))
                fm = getRegression(tmtA.cov.vars.df,"tmtA","glm","ALL")
            
            }else if (run_type == 'tmtB'){
                Y.tmtB.vec = func.tmtB.cgt(cgtData,mapData,varMat)
                Y.tmtB.log = log(Y.tmtB.vec)
                tmtB.cov.vars.df = data.frame(cbind(Y.tmtB.log,varMat.cov[,!(colnames(varMat.cov) %in% "f.eid")]))
                fm = getRegression(tmtB.cov.vars.df,"tmtB","glm","ALL")
            }else if (run_type == 'fit'){
                Y.fit.vec = func.fit.cgt(cgtData,mapData,varMat)
                fit.cov.vars.df = data.frame(cbind(Y.fit.vec,varMat.cov[,!(colnames(varMat.cov) %in% "f.eid")]))
                fm = getRegression(fit.cov.vars.df,"fit","glm","ALL")
            
            }else if (run_type == 'dst'){
                Y.dst.vec = func.dst.cgt(cgtData,mapData,varMat)
                dst.cov.vars.df = data.frame(cbind(Y.dst.vec,varMat.cov[,!(colnames(varMat.cov) %in% "f.eid")]))
                fm = getRegression(dst.cov.vars.df,"dst","glm","ALL")

            }else if (run_type == 'pmt'){
                Y.pmt.vec = func.pmt.cgt(cgtData,mapData,varMat)
                Y.pmt.log = log(1+Y.pmt.vec)
                pmt.cov.vars.df = data.frame(cbind(Y.pmt.log,varMat.cov[,!(colnames(varMat.cov) %in% "f.eid")]))
                fm = getRegression(pmt.cov.vars.df,"pmt","glm","ALL")
            }
            
            save(fm,file=fm_obj_path)

            res.df = data.frame(summary(fm)$coefficients)
            res.df[,4] = as.numeric(res.df[,4])
            res.df.sort = res.df[order(res.df[,4],decreasing=F),]
            res.df.sort.2 = data.frame(cbind(rownames(res.df.sort),res.df.sort))
            colnames(res.df.sort.2) = c("Variable","Estimate","StdError","T-Value","P-Value")
            write.table(res.df.sort.2,sep='\t',file=sort_sig_path,col.names=T,row.names=F,quote=F)
        
        }else if(run_type %in% c("qualf","occup","hhi")) {
            
            X_alc_df = getAlcoDF(cgtData,mapData,varMat)
            #varMat.cov = data.frame(cbind(X_sex_df,X_age_df,X_alc_df,data.frame(varMat)))
            varMat.cov = data.frame(cbind(X_sex_df,X_age_df,data.frame(varMat)))

            if (run_type == 'qualf'){
                Y.qual.vec = getQualDF(gmfData,mapData,varMat)
                reg.raw.df.tmp = cbind(Y.qual.vec,varMat.cov)
                reg.raw.df = reg.raw.df.tmp[,-which(colnames(reg.raw.df.tmp)=='f.eid')]
                colnames(reg.raw.df)[1] = 'Qualifications'

                #Filter the regression DF for Qualification code= -3 and -Inf
                reg.filt.1 = reg.raw.df[reg.raw.df$Qualifications!=-3,]
                reg.filt.2 = reg.filt.1[!is.infinite(reg.filt.1$Qualifications),]
                reg.filt.df = reg.filt.2

                #Assign Qualification code
                reg.fc.df = data.frame(getQualCode(reg.filt.df))
                #reg.fc.df$Qualifications = factor(reg.fc.df$Qualifications)
                reg.fc.df$Qualifications = ordered(reg.fc.df$Qualifications,
                                                   levels=c("A-Level","College","CSE",
                                                            "NVQ","O-Level","Other","NOTA"
                                                            )
                                                   )

            }else if(run_type =="occup"){
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
           
            stop()
            fm = getRegression(reg.fc.df,run_type,'ordinal','ALL')

            ctable <- coef(summary(fm))
            p = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
            ctable <- cbind(ctable, "p value" = p)

            #ci <- confint(fm)
            ci.default = confint.default(fm)
            lor = exp(coef(fm))
            lor_ci = data.frame(exp(cbind(OR = coef(fm), ci.default)))
            #lor_ci.df = cbind(rownames(lor_ci),lor_ci)

            or_vec = rep('NA',dim(ctable)[1])
            or_25 = rep('NA',dim(ctable)[1])
            or_95 = rep('NA',dim(ctable)[1])

            ctable.lor.df = data.frame(cbind(ctable,or_vec,or_25,or_95))
            rownames(ctable.lor.df) = rownames(ctable)
            

            for (i in 1:dim(ctable)[1]){
                r_name = rownames(ctable)[i]
                ind_lor = which(rownames(lor_ci)==r_name)
                if(length(ind_lor !=0)) {
                    ctable.lor.df[i,5:7] = as.vector(lor_ci[ind_lor,] )
                }
            }
            
            res.df = ctable.lor.df
            res.df$p.value = as.numeric(res.df$p.value)
            res.df.sort = res.df[order(res.df[,'p.value'],decreasing=F),]
            res.df.sort.2 = cbind(rownames(res.df.sort),res.df.sort)

            colnames(res.df.sort.2) = c("Predictors","Estimate","StdError","T-Value","P-Value","OR","2.5%","97.5%")
            write.table(res.df.sort.2,sep='\t',file=sort_sig_path,col.names=T,row.names=F,quote=F)

            result = list()
            result[[1]] = reg.fc.df
            result[[2]] = fm
            result[[3]] = ctable
            result[[4]] = lor_ci
            result[[5]] = res.df

            save(result,file=fm_obj_path)
            ## Printing LOR output
            #colnames(lor_ci.df) = c("Variables","OR","2.5%","97.5%")
            #sort_lor_path = paste0(out_res_dir,"/",gene_id,'_sort_lor.txt',collapse='')
            #write.table(lor_ci.df,sep='\t',file=sort_lor_path,col.names=T,row.names=F,quote=F)
        }
    }
}
