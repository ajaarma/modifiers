args=(commandArgs(TRUE))
loc_4 = "/hps/software/users/dunham/R_lib/4.0.3/"
require('data.table',lib=loc_4)
require('Matrix',lib=loc_4)
require('moments',lib=loc_4)

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

batchNum = 'merge_all'
#test_str = c("cgt_pmt_res","cgt_fit_res","cgt_dst_res","cgt_sds_res","cgt_tmtA_res","cgt_tmtB_res")

run_type_list = args[[1]]
file_str = args[[2]]

id_int_type = 'int_genes_results'
id_int_merge_type = 'int_genes_ukbb_merge_all'

#run_type_list = c('pmt','fit','dst','sds','tmtA','tmtB')

for (run_type in run_type_list) {
    if (run_type =='sds') {
        y_ind = c(73)
        test_str = c("cgt_sds_res")
    }else if(run_type =='tmtA'){
        y_ind = c(71)
        test_str = c("cgt_tmtA_res")
    }else if(run_type =='tmtB'){
        y_ind = c(72)
        test_str = c("cgt_tmtB_res")
    }else if(run_type =='fit'){
        y_ind = c(74)
        test_str = c("cgt_fit_res")
    }else if(run_type =='dst'){
        y_ind = c(75)
        test_str = c('cgt_dst_res')
    }else if(run_type =='pmt'){
        y_ind = c(68:70)
        test_str = c('cgt_pmt_res')
    }

    for(test_id in test_str) {
        
        out_res_dir = paste("/nfs/research/dunham/samples/ukbb/data/results",id_int_type,test_id,sep="/")
        system(paste0('mkdir -p ',out_res_dir,collapse=""))
        
        log_file = paste0(out_res_dir,'/.',file_str,'_log_file.txt',collapse='')
        sink(log_file)
        
        #file_str = c("_vt_rc25_hc6","_vt_urcgte25_hc6","_vt_urclte25_hc6")
        #file_str = c("_vt_urcgte25_hc6","_vt_urclte25_hc6")

        for (file_id in file_str) {
            
            #ukb_var_path = paste0("/hps/nobackup/dunham/ukbb/filterSNV/ukbb_",batchNum,"/id_genes/variantsTable.txt",collapse="")
            #ukb_var_path = paste0("/hps/nobackup/dunham/ukbb/filterSNV/ukbb_",batchNum,"/id_genes/variantsTable",file_id,".txt",collapse="")
            #ukb_var_path = paste0("/hps/nobackup/dunham/ukbb/filterSNV/ukbb_",batchNum,'/merge_all',file_id,".txt",collapse="")
            ukb_var_path = paste0("/hps/nobackup/dunham/ukbb/filterSNV/",id_int_merge_type,'/merge_all',file_id,".txt",collapse="")
            ukb_map_path = paste0("/nfs/research/dunham/samples/ukbb/data/icd10/",batchNum,"_map.txt",collapse="")
            ukb_ex_sample = "/nfs/research/dunham/samples/ukbb/data/icd10/NonNeuroF7_samples_50K.txt"
            ukb_test_path = "/nfs/research/dunham/samples/ukbb/data/cg_test/samples_cg_test_v3.txt"
            sigRes_path = paste0("/nfs/research/dunham/samples/ukbb/data/results/",id_int_type,"/",test_id,
                                 "/",batchNum,file_id,'_sig.txt',collapse='')
            sigResSort_path = paste0("/nfs/research/dunham/samples/ukbb/data/results/",id_int_type,"/",test_id,
                                     "/",batchNum,file_id,'_sort_sig.txt',collapse='')

            #Read UKBB non-ID
            ukbNonID = fread(ukb_ex_sample,sep='\t',header=F,stringsAsFactor=F,quote="")
            id_samples = ukbNonID$V2

            #Read UKB variants
            varData = fread(ukb_var_path,sep="\t",header=F,stringsAsFactor=F,quote="")

            #Read UKB Ma file
            mapData = fread(ukb_map_path,sep="\t",header=F,stringsAsFactor=F,quote="")

            #Read UKB Cognitive Test
            cgtData = fread(ukb_test_path,sep="\t",header=T,stringsAsFactor=F,quote="")

            varMat =  getVarMat(varData,mapData,id_samples)
           
            # Pairs-Matching Test
            #x_ind = c(1,2,76:79);y_ind = c(68:70)
            x_ind = c(1,2,76:79);
            
            pairs_test_list = getCgtTest(cgtData,mapData,rownames(varMat),x_ind, y_ind)#(1:6,31:33))
            map_sort_df = pairs_test_list[[1]]
            cgt_X_pmt_df = pairs_test_list[[2]]

            # Pairs-Matching Test
            cgt_Y_pmt_df = getPMTvector(cgt_X_pmt_df,map_sort_df,run_type)
            
            #Regression analysis
            gfm_list = list()
            geneList = as.vector(unique(factor(varData$V5)))
            
            sigRes = data.frame()

            for (i in 1:length(geneList)){
                gene_id = geneList[i]
                
                X_mat = cbind(varMat,cgt_X_pmt_df[,c(2:4)])
                cat(i,'\t',gene_id,'\t',file_id,'\t',run_type,'\n')
                
                gfm_list[[gene_id]] = getGLM(varMat,cgt_X_pmt_df,cgt_Y_pmt_df,gene_id,run_type)
                a11 = data.frame(summary(gfm_list[[gene_id]])$coefficients)
                a11.rows = gsub('^X','',rownames(a11))
                row_prefix = rep(gene_id,length(a11.rows))
                a11_prefix.rows = paste0(row_prefix,'_',a11.rows,sep="")

                #tmpSigRes = data.frame(coef(summary(gfm_list[[gene_id]]))[,4])
                #sigRes = rbind(sigRes,tmpSigRes)
                a11.tmp = cbind(a11_prefix.rows,a11)
                sigRes = rbind(sigRes,a11.tmp)

                gfm_list = list()
            }
            
            sigRes.uniq = unique(sigRes)
            sigRes.uniq.vars = cbind(rownames(sigRes.uniq),sigRes.uniq)
            colnames(sigRes.uniq) = c("Variables","Estimate","Std.Error","T-value","FDRsignificance")
            sigRes.uniq.sort = sigRes.uniq[order(sigRes.uniq$FDRsignificance),]
            write.table(sigRes.uniq.sort,sep='\t',file=sigResSort_path,row.names=F,col.names=T,quote=F)
            write.table(sigRes.uniq,sep='\t',file=sigRes_path,row.names=F,col.names=T,quote=F)
        }
    }
}
# Generalized Linear model

#c("SLC26A4" "SBDS"    "PEX1"    "GLI3"    "FAM20C"  "ASL"     "UPF1"
# [8] "PNKP"    "NDUFS7"  "GAMT"    "DOCK6"   "CERS1"   "USH2A"   "TSHB"
#[15] "TRIT1"   "POMGNT1" "NEK2"    "HPDL"    "HAX1"    "DPYD"    "CRB1"
#[22] "ABCA4"


#BAyesian Linear regression
#DF <-5
#Vy = var(cg_test_Y$PMT_Score,na.rm=T)
#h2 = 0.5
#Se = Vy*(1-h2)*(DF-2)
#MSx = sum(apply(FUN=var,MARGIN=2,X=varMat))
#Sr = Vy*h2*(DF-2)/MSx
#prior = list(varE=list(df=DF,S=Se),varBR=list(df=DF,S=Sr))

#fm = BLR(y=cg_test_Y$PMT_Score,XF=varMat,nIter=20500,burnIn=2000,prior=prior)
#lm.fm = lm(cg_test_Y$PMT_Score~X7_107674970_SLC26A4+X7_107675051_SLC26A4+X7_107661637_SLC26A4+X7_107689054_SLC26A4,data=varMat)
