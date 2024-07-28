args=(commandArgs(TRUE))
if(Sys.info()['sysname']!='Darwin') {
    loc_4 = "/hps/software/users/dunham/R_lib/4.0.3/"
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
}else{
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
  require('ggplot2')
}

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


getVarMat = function(df,mapData,id_samples){

    df = data.frame(df)
    ukb_samples = data.frame(mapData)[,1]

    #row_names = unique(as.vector(df[,2]))
    row_names = ukb_samples[!ukb_samples %in% id_samples]
    col_names = unique(as.vector(df[,1]))

    mat = mat.or.vec(length(row_names),length(col_names))
    colnames(mat) = col_names
    rownames(mat) = row_names

    for (i in 1:dim(df)[1]){
        row_id = df[i,2]
        col_id = df[i,1]
        r_index = which(rownames(mat)==row_id)
        c_index = which(colnames(mat)==col_id)
        mat[r_index,c_index] = 1
    }
    return(mat)
}

getGender <- function(gmfData,mapData,ukb_var_col,gender_field){
    
    #Subroutine to extract covariates for non-ID samples
    
    y_ind = c(which(grepl('^f.31.',colnames(gmfData))==TRUE))

    # Map and sort according to UKB sample ids
    ukb_eid_df = subset(mapData,mapData$V1 %in% ukb_var_col)
    ukb_eid_df2 = ukb_eid_df[match(ukb_var_col,ukb_eid_df$V1),]

    #Gender
    sex_df = data.frame(gmfData)[,y_ind]
    #sex.rfm.df = reformGender(sex_df)
    sex.rfm.df = data.frame(sex_df)
    
    sex.df = data.frame(cbind(gmfData$f.eid,sex.rfm.df))
    colnames(sex.df) = c('f.eid','Gender')
    #colnames(sex.df) = c('f.eid')

    X_sex_df2 = subset(sex.df,sex.df$f.eid %in% ukb_eid_df2$V2)
    X_sex_df = X_sex_df2[match(ukb_eid_df2$V2,X_sex_df2$f.eid),]

    return(X_sex_df)

}

getQualDF <- function(gmfData,mapData,ukb_var_col,field_id_ind){

    #### Get Independent variable Qualification DF ###

    y_ind = c(which(grepl('^f.6138.',colnames(gmfData))==TRUE))

    # Map and sort according to UKB sample ids
    ukb_eid_df = subset(mapData,mapData$V1 %in% ukb_var_col)
    ukb_eid_df2 = ukb_eid_df[match(ukb_var_col,ukb_eid_df$V1),]

    #Qualifications
    qual_df = data.frame(gmfData)[,y_ind]
    
    qual.df = cbind(gmfData$f.eid,qual_df)
    colnames(qual.df)[1] = c('f.eid')

    Y_qual_df2 = subset(qual.df,qual.df$f.eid %in% ukb_eid_df2$V2)
    Y_qual_df = Y_qual_df2[match(ukb_eid_df2$V2,Y_qual_df2$f.eid),]

    Y.qual.vec = apply(Y_qual_df,1,function(x){
                                tmp = min(x[-1],na.rm=T)
                                return(tmp)
                        })

    names(Y.qual.vec) = as.vector(Y_qual_df[,1])
    #colnames(Y.qual.vec) = c("QualF")

    return(Y.qual.vec)
}


getQualCode <- function(reg.filt.df){

    Y.qual.tmp = reg.filt.df
    Y.qual.tmp[Y.qual.tmp$Qualifications==1,1] = 'College'
    Y.qual.tmp[Y.qual.tmp$Qualifications==2,1] = 'A-Level'
    Y.qual.tmp[Y.qual.tmp$Qualifications==3,1] = 'O-Level'
    Y.qual.tmp[Y.qual.tmp$Qualifications==4,1] = 'CSE'
    Y.qual.tmp[Y.qual.tmp$Qualifications==5,1] = 'NVQ'
    Y.qual.tmp[Y.qual.tmp$Qualifications==6,1] = 'Other'
    Y.qual.tmp[Y.qual.tmp$Qualifications==-7,1] = 'NOTA'

    return(Y.qual.tmp)
}
getOccupCode <- function(reg.filt.df){

    Y.occup.tmp = reg.filt.df
    Y.occup.tmp[,1] = apply(data.frame(Y.occup.tmp[,1]),1,function(x){
                            tmp = as.numeric(strsplit(as.character(x),"")[[1]][1])
                            return(tmp)
                        })
    colnames(Y.occup.tmp)[1] = 'Occupation'
    Y.occup.tmp[Y.occup.tmp$Occupation==1,1] = 'SRManager'
    Y.occup.tmp[Y.occup.tmp$Occupation==2,1] = 'ProfOccup'
    Y.occup.tmp[Y.occup.tmp$Occupation==3,1] = 'AssocProfTechOccup'
    Y.occup.tmp[Y.occup.tmp$Occupation==4,1] = 'AdminOccup'
    Y.occup.tmp[Y.occup.tmp$Occupation==5,1] = 'SkillTradeOccup'
    Y.occup.tmp[Y.occup.tmp$Occupation==6,1] = 'SelfServiceOccup'
    Y.occup.tmp[Y.occup.tmp$Occupation==7,1] = 'SalesCustomOccup'
    Y.occup.tmp[Y.occup.tmp$Occupation==8,1] = 'ProcPlantMachineOpt'
    Y.occup.tmp[Y.occup.tmp$Occupation==9,1] = 'ElementaryOccup'

    return(Y.occup.tmp)
}

getHHICode <- function(reg.filt.df){

    Y.hhi.tmp = reg.filt.df[!is.na(reg.filt.df$HHI),]
    
    Y.hhi.tmp[Y.hhi.tmp$HHI==1,1] = 'LT-18K'
    Y.hhi.tmp[Y.hhi.tmp$HHI==2,1] = '18K-31K'
    Y.hhi.tmp[Y.hhi.tmp$HHI==3,1] = '31K-52K'
    Y.hhi.tmp[Y.hhi.tmp$HHI==4,1] = '52K-100K'
    Y.hhi.tmp[Y.hhi.tmp$HHI==5,1] = 'GT-100K'
    
    return(Y.hhi.tmp)
}


getOccupDF <- function(gmfData,mapData,ukb_var_col,field_id_ind){

    #### Get Independent variable Occupation DF ###

    y_ind = c(which(grepl('^f.22617.',colnames(gmfData))==TRUE))

    # Map and sort according to UKB sample ids
    ukb_eid_df = subset(mapData,mapData$V1 %in% ukb_var_col)
    ukb_eid_df2 = ukb_eid_df[match(ukb_var_col,ukb_eid_df$V1),]

    #Qualifications
    occup_df = data.frame(gmfData)[,y_ind]
    
    occup.df = cbind(gmfData$f.eid,occup_df)
    colnames(occup.df)[1] = c('f.eid')

    Y_occup_df2 = subset(occup.df,occup.df$f.eid %in% ukb_eid_df2$V2)
    Y_occup_df = Y_occup_df2[match(ukb_eid_df2$V2,Y_occup_df2$f.eid),]

    Y.occup.vec = apply(Y_occup_df,1,function(x){
                                x = x[-1]
                                tmp = x[!is.na(x)]
                                return(min(tmp))
                            })

    names(Y.occup.vec) = as.vector(Y_occup_df[,1])
    #colnames(Y.qual.vec) = c("QualF")

    return(Y.occup.vec)
}
getHhiDF <- function(gmfData,mapData,ukb_var_col,field_id_ind){

    #### Get Independent variable House Hold Income (HHI) DF ###

    y_ind = c(which(grepl('^f.738.',colnames(gmfData))==TRUE))

    # Map and sort according to UKB sample ids
    ukb_eid_df = subset(mapData,mapData$V1 %in% ukb_var_col)
    ukb_eid_df2 = ukb_eid_df[match(ukb_var_col,ukb_eid_df$V1),]

    #HHI
    hhi_df = data.frame(gmfData)[,y_ind]
    
    hhi.df = cbind(gmfData$f.eid,hhi_df)
    colnames(hhi.df)[1] = c('f.eid')

    Y_hhi_df2 = subset(hhi.df,hhi.df$f.eid %in% ukb_eid_df2$V2)
    Y_hhi_df = Y_hhi_df2[match(ukb_eid_df2$V2,Y_hhi_df2$f.eid),]
    Y.hhi.vec = Y_hhi_df[,2]

    names(Y.hhi.vec) = as.vector(Y_hhi_df[,1])

    return(Y.hhi.vec)
}


getRegDF <- function(reg.fc.df,gene_id,rest_cov_index,run_type) {
    ##### Function to get Ordinal Regression Data frame #####
    
    gene_id.grep = paste("_",gene_id,"$",sep="")
    #gene_id.grep = paste("^",gene_id,"_",sep="")
    index = which(grepl(gene_id.grep,colnames(reg.fc.df))==T)

    if (run_type == 'qualf'){
        reg.gene.run = as.data.frame(reg.fc.df$Qualifications)
        colnames(reg.gene.run) = "Qualifications"
        reg.gene.run$Qualifications = factor(reg.gene.run$Qualifications)
    }else if(run_type =="occup"){
        reg.gene.run = as.data.frame(reg.fc.df$Occupation)
        colnames(reg.gene.run) = "Occupation"
        reg.gene.run$Occupation = factor(reg.gene.run$Occupation)
    }else if(run_type == 'hhi'){
        reg.gene.run = as.data.frame(reg.fc.df$HHI)
        colnames(reg.gene.run) = "HHI"
        reg.gene.run$HHI = factor(reg.gene.run$HHI)
    }
    
    reg.gene.vars = as.data.frame(reg.fc.df[,index])
    colnames(reg.gene.vars) = colnames(reg.fc.df)[index]

    reg.gene.cov = as.data.frame(reg.fc.df[,rest_cov_index])
    colnames(reg.gene.cov) = colnames(reg.fc.df)[rest_cov_index]

    #reg.gene.df = data.frame(cbind(reg.fc.df$Qualifications,reg.fc.df[,index],reg.fc.df[,rest_cov_index]))
    reg.gene.df = data.frame(cbind(reg.gene.run,reg.gene.vars,reg.gene.cov))
    
    #colnames(reg.gene.df)[1] = "Qualifications"
    #colnames(reg.gene.df)[dim(reg.gene.df)[2]] = "Gender"
    #reg.gene.df$Qualifications = factor(reg.gene.df$Qualifications)
    
    fmt_str_lhs = colnames(reg.gene.df)[1]
    fmt_str_rhs = paste0(colnames(reg.gene.df)[-1],collapse="+")

    formula_str = paste0(fmt_str_lhs," ~ ",fmt_str_rhs,collapse="")

    fm <- polr(formula_str,data=reg.gene.df,Hess=TRUE)
    
    return(fm)

}

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
#test_str = c("cgt_pmt_res","cgt_fit_res","cgt_dst_res","cgt_sds_res","cgt_tmtA_res","cgt_tmtB_res")

run_type_list = 'qualf' #args[[1]]
file_str =  '_vt_rc25_hc6' #args[[2]]
id_int_str =  'id' #args[[3]]
gene_id =  'SLC26A4' #args[[4]]
out_res_dir = '/hps/nobackup/dunham/ukbb/id_genes/occup/_vt_rc25_hc6' #args[[5]]

if (id_int_str=="int"){
    id_int_type = 'int_genes_results'
    id_int_merge_type = 'int_genes_ukbb_merge_all'
}else if (id_int_str == 'id'){
    id_int_type = 'id_genes_results'
    id_int_merge_type = 'id_genes_ukbb_merge_all'
}

#UKB GMF data
ukb_gmf_path = paste0(ukb_path,"/data/gmf_test/samples_gmf_test.txt",collapse="")
gmfData = fread(ukb_gmf_path,sep="\t",header=T,stringsAsFactor=F,quote="")
#gmfData = data.frame(gmfData)
gmf_col = colnames(gmfData)

#UKB ICD10 Map path
ukb_map_path = paste0(ukb_path,"/data/icd10/",batchNum,"_map.txt",collapse="")
mapData = fread(ukb_map_path,sep="\t",header=F,stringsAsFactor=F,quote="")
#mapData = data.frame(mapData)

#UKB non-ID samples
ukb_ex_sample = paste0(ukb_path,"/data/icd10/NonNeuroF7_samples_50K.txt",collapse="")
ukbNonID = fread(ukb_ex_sample,sep='\t',header=F,stringsAsFactor=F,quote="")
#ukbNonId = data.frame(ukbNonID)
id_samples = ukbNonID$V2


#Extract col index per GMF type
for (run_type in run_type_list) {
    
    #Qualification
    if (run_type =='qualf') {
        y_ind =  which(grepl('f.6138.',gmf_col)==TRUE)
        test_str = c("gmf_qualf_res")
        gmf.df = getGMF(gmfData,run_type,y_ind)

    #Occupation
    }else if(run_type =='occup'){
        y_ind = which(grepl('f.22617.',gmf_col)==TRUE)
        test_str = c("gmf_occup_res")

    #Houseghold Income
    }else if(run_type =='hhi'){
        y_ind = which(grepl('f.738.',gmf_col)==TRUE)
        test_str = c("gmf_hhi_res")

    #Townsend Deprivation Index
    }else if(run_type =='tdi'){
        y_ind = which(grepl('f.457.',gmf_col)==TRUE)
        tdi_x_col = which(grepl('f.680|f.728|f.709|f.6142',gmf_col)==TRUE)
        test_str = c("gmf_tdi_res")
    }

    for(test_id in test_str) {
        
        for (file_id in file_str) {
            
            ukb_var_path = paste0(ukb_path,"/analysis/filterSNV/",id_int_merge_type,'/merge_all',file_id,".txt",collapse="")
            ctable_path = paste0(out_res_dir,"/",gene_id,file_id,'_ctable.txt',collapse='')
            lor_path = paste0(out_res_dir,"/",gene_id,file_id,'_lor.txt',collapse='')

            #Read UKB variants
            varData = fread(ukb_var_path,sep="\t",header=F,stringsAsFactor=F,quote="")
            
            #Get variant-sample Matrix (rows=UKBB non-ID, col=variants)
            varMat =  getVarMat(varData,mapData,id_samples)
            #col.name.vec = apply(data.frame(colnames(varMat)),1,function(x){
            #                                strs = strsplit(x,"_")[[1]]
            #                                col_str = paste0(strs[3],"_",
            #                                                 strs[1],"_",
            #                                                 strs[2],collapse="")
            #                                return(col_str)
            #                                }
            #                    )
            #colnames(varMat) = col.name.vec

            #### Covariates to be included ###############
            # Gender X-variable
            gender_field = 'f.31'

            #Rent Acco, Num of House-hold, Num of Cars, Current emp status
            rent_acco_field = 'f.680'
            hh_num_field = 'f.709'
            car_num_field = 'f.728'
            emp_stat_field = 'f.6142'

            X_sex_df = getGender(gmfData,mapData,rownames(varMat),gender_field) 

            #Returns vector with scores: -7,-3,1,2,3,4,5,6,-Inf
            #-Inf: Corresponds to those individuals with all NAs for each
            #instances
            if (run_type == 'qualf') {
                Y_qual_vec = getQualDF(gmfData,mapData,rownames(varMat),gender_field)
                reg.raw.df.tmp = cbind(Y_qual_vec,varMat,X_sex_df)
                reg.raw.df = reg.raw.df.tmp[,-which(colnames(reg.raw.df.tmp)=='f.eid')]
                colnames(reg.raw.df)[1] = 'Qualifications'

                #Filter the regression DF for Qualification code= -3 and -Inf
                reg.filt.1 = reg.raw.df[reg.raw.df$Qualifications!=-3,]
                reg.filt.2 = reg.filt.1[!is.infinite(reg.filt.1$Qualifications),]
                reg.filt.df = reg.filt.2

                #Assign Qualification code
                reg.fc.df = data.frame(getQualCode(reg.filt.df))
                #rest_cov_index = c(which(colnames(reg.fc.df)=='Female'), which(colnames(reg.fc.df)=='Male'))
                rest_cov_index = c(which(colnames(reg.fc.df)=='Gender'))

            }else if(run_type=='occup'){
                Y_occup_vec = getOccupDF(gmfData,mapData,rownames(varMat),gender_field)
                reg.raw.df.tmp = cbind(Y_occup_vec,varMat,X_sex_df)
                reg.raw.df = reg.raw.df.tmp[,-which(colnames(reg.raw.df.tmp)=='f.eid')]
                colnames(reg.raw.df)[1] = 'Occupation'

                #Filter the regression Occup-DF for Occupation -Inf and apply
                #Occuptation codes
                reg.filt.1 = reg.raw.df[!is.infinite(reg.raw.df$Occupation),]
                reg.filt.2 = data.frame(getOccupCode(reg.filt.1))
                reg.fc.df = reg.filt.2
                rest_cov_index = c(which(colnames(reg.fc.df)=='Gender'))
            
            }else if (run_type=='hhi'){
                Y_hhi_vec = getHhiDF(gmfData,mapData,rownames(varMat),gender_field)
                reg.raw.df.tmp = cbind(Y_hhi_vec,varMat,X_sex_df)
                reg.raw.df = reg.raw.df.tmp[,-which(colnames(reg.raw.df.tmp)=='f.eid')]
                colnames(reg.raw.df)[1] = 'HHI'

                #Filter the regression HHI-DF for Occupation -Inf and apply
                #Occuptation codes
                reg.filt.1 = reg.raw.df[reg.raw.df$HHI !=-1,]
                reg.filt.2 = reg.filt.1[reg.filt.1$HHI !=-3,]
                reg.fc.df = data.frame(getHHICode(reg.filt.2))
                rest_cov_index = c(which(colnames(reg.fc.df)=='Gender'))
            
            }
            ct.df  =  data.frame()
            lor.df = data.frame()

            m = getRegDF(reg.fc.df,gene_id,rest_cov_index,run_type)
            stop() 
            m_obj_file = paste0(out_res_dir,'/',gene_id,'_reg_obj.RData',collapse="")
            save(m,file=m_obj_file)

            ctable <- coef(summary(m))
            p = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
            ctable <- cbind(ctable, "p value" = p)
            ci <- confint(m)
            ci.default = confint.default(m)
            lor = exp(coef(m))
            lor_ci = data.frame(exp(cbind(OR = coef(m), ci)))

            a11.rows = rownames(ctable)
            row_prefix = rep(gene_id,dim(ctable)[1])
            ct_row_prefix = paste0(row_prefix,'_',a11.rows,sep="")
            ct.tmp = cbind(ct_row_prefix,ctable)

            b11.rows = rownames(lor_ci)
            row_prefix = rep(gene_id,dim(lor_ci)[1])
            lor_row_prefix = paste0(row_prefix,'_',b11.rows,sep="")
            lor.tmp = cbind(lor_row_prefix,lor_ci)

            ct.df = rbind(ct.df,ct.tmp)
            lor.df = rbind(lor.df,lor.tmp)
            colnames(ct.df)= c("Predictors","Value","Std.error","T-value","P-value")
            #cat(gene_id,'\n',file=sigRes_path)
            write.table(ct.df,sep='\t',file=ctable_path,row.names=F,col.names=T,quote=F)

            #Printing LOR output
            colnames(lor.df)= c("Predictors","OR","2.5%","97.5%")
            write.table(lor.df,sep='\t',file=lor_path,row.names=F,col.names=T,quote=F)

        }
    }
}
