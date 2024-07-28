args=(commandArgs(TRUE))
if(Sys.info()['sysname']!='Darwin') {
    loc_4 = "/hps/software/users/dunham/R_lib/4.0.3/"
    require('R.utils')
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
    require('ltm')
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
  require('ltm')
}

sourceDirectory(paste("ukbb/",sep=""))
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

ukb_path = "/nfs/research/dunham/samples/ukbb/"
batchNum = 'merge_all'

####### General Variable ###############

run_type_list = 'tdi' #args[[1]]
file_str =  '_vt_rc25_hc6' #args[[2]]
out_res_dir = '/hps/nobackup/dunham/ukbb/id_genes/occup/_vt_rc25_hc6' #args[[5]]


#UKB CGT data
ukb_cgt_path = paste0(ukb_path,"/data/cg_test/samples_cg_test_v3.txt",collapse="")
cgtData = fread(ukb_cgt_path,sep="\t",header=T,stringsAsFactor=F,quote="")
cgt_col = colnames(cgtData)
message('Finished reading CGT data')



#UKB GMF data
ukb_gmf_path = paste0(ukb_path,"/data/gmf_test/samples_gmf_test.txt",collapse="")
gmfData = fread(ukb_gmf_path,sep="\t",header=T,stringsAsFactor=F,quote="")
#gmfData = data.frame(gmfData)
gmf_col = colnames(gmfData)
message('Finished reading GMF data')

#UKB ICD10 Map path
ukb_map_path = paste0(ukb_path,"/data/icd10/",batchNum,"_map.txt",collapse="")
mapData = fread(ukb_map_path,sep="\t",header=F,stringsAsFactor=F,quote="")
#mapData = data.frame(mapData)
message('Finished reading UKB-ICD10 map data')

#UKB non-ID samples
ukb_ex_sample = paste0(ukb_path,"/data/icd10/NonNeuroF7_samples_50K.txt",collapse="")
ukbNonID = fread(ukb_ex_sample,sep='\t',header=F,stringsAsFactor=F,quote="")
#ukbNonId = data.frame(ukbNonID)
id_samples = ukbNonID$V2
message('Finished reading UKB-nonID samples')

#Extract col index per GMF type
#ukb_var_path = paste0(ukb_path,"/analysis/filterSNV/",id_int_merge_type,'/merge_all',file_id,".txt",collapse="")
ukb_var_path = "/nfs/research/dunham/samples/ukbb/analysis/filterSNV/id_genes_ukbb_merge_all/merge_all_vt_rc25_hc6.txt"

#Read UKB variants
varData = fread(ukb_var_path,sep="\t",header=F,stringsAsFactor=F,quote="")
message('Finished reading UKB variants list')

#Get variant-sample Matrix (rows=UKBB non-ID, col=variants)
varMat =  getVarMat(varData,mapData,id_samples)

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

for (run_type in run_type_list) {
    if (run_type == 'qualf') {
        Y_qual_vec = getQualDF(gmfData,mapData,rownames(varMat),'qualf')

    }else if(run_type=='occup'){
        Y_occup_vec = getOccupDF(gmfData,mapData,rownames(varMat),'occup')

    }else if (run_type=='hhi'){
        Y_hhi_vec = getHhiDF(gmfData,mapData,rownames(varMat),'hhi')

    }else if (run_type =='tdi'){

        Y_tdi_vec = getTdiDF(gmfData,mapData,varMat,'tdi')
        X_ncars_vec = getTdiDF(gmfData,mapData,varMat,'ncars')
        X_nhh_vec = getTdiDF(gmfData,mapData,varMat,'nhh')
        X_rent_vec = getTdiDF(gmfData,mapData,varMat,'rent')
        X_emp_vec = getTdiDF(gmfData,mapData,varMat,'emp')

    }
}



tdi.comb.df = as.data.frame(cbind(Y_tdi_vec,X_ncars_vec,X_nhh_vec,X_rent_vec,X_emp_vec))
tdi.comb.df.2 = tdi.comb.df[!is.na(tdi.comb.df$tdi),]

tdi.ncars.outList = func.tdi.ncars.cor(tdi.comb.df.2)
tdi.ncars.df = tdi.ncars.outList[[1]]
tdi.ncars.cor = tdi.ncars.outList[[2]]

tdi.nhh.outList = func.tdi.nhh.cor(tdi.comb.df.2)
tdi.nhh.df = tdi.nhh.outList[[1]]
tdi.nhh.cor = tdi.nhh.outList[[2]]

tdi.rent.outList = func.tdi.rent.cor(tdi.comb.df.2)
tdi.rent.df = tdi.rent.outList[[1]]
tdi.rent.cor = tdi.rent.outList[[2]]

tdi.emp.outList = func.tdi.emp.cor(tdi.comb.df.2)
tdi.emp.df = tdi.emp.outList[[1]]
tdi.emp.cor = tdi.emp.outList[[2]]

########### Correlation between Cognitive-Tests and TDI ###############
Y.sds.vec = func.sds.cgt(cgtData,mapData,varMat)
tdi.sds.cor = cor(tdi.comb.df$tdi,Y.sds.vec,use="na.or.complete")


Y.tmtA.vec = func.tmtA.cgt(cgtData,mapData,varMat)
(tdi.tmtA.cor = cor(tdi.comb.df$tdi,Y.tmtA.vec,use="na.or.complete"))

Y.tmtB.vec = func.tmtB.cgt(cgtData,mapData,varMat)
(tdi.tmtB.cor = cor(tdi.comb.df$tdi,Y.tmtB.vec,use="na.or.complete"))

Y.fit.vec = func.fit.cgt(cgtData,mapData,varMat)
(tdi.fit.cor = cor(tdi.comb.df$tdi,Y.fit.vec,use="na.or.complete"))

Y.dst.vec = func.dst.cgt(cgtData,mapData,varMat)
(tdi.dst.cor = cor(tdi.comb.df$tdi,Y.dst.vec,use="na.or.complete"))

Y.pmt.vec = func.pmt.cgt(cgtData,mapData,varMat)
(tdi.pmt.cor = cor(tdi.comb.df$tdi,Y.pmt.vec,use="na.or.complete"))


########### Correlation Matrix ##################
cgt.gmf.common = Reduce(intersect,list(rownames(Y_tdi_vec), rownames(tdi.ncars.df),
                                       rownames(tdi.nhh.df), rownames(tdi.rent.df),
                                       rownames(tdi.emp.df),rownames(Y.sds.vec),
                                       rownames(Y.tmtA.vec),rownames(Y.tmtB.vec),rownames(Y.fit.vec),
                                       rownames(Y.dst.vec),rownames(Y.pmt.vec)
                                   )
                    )

Y.tdi.vec.df2 = Y_tdi_vec[match(cgt.gmf.common,rownames(Y_tdi_vec)),]
tdi.ncars.df2 = tdi.ncars.df[match(cgt.gmf.common,rownames(tdi.ncars.df)),]
tdi.nhh.df2 = tdi.nhh.df[match(cgt.gmf.common,rownames(tdi.nhh.df)),]
tdi.rent.df2 = tdi.rent.df[match(cgt.gmf.common,rownames(tdi.rent.df)),]
tdi.emp.df2 = tdi.emp.df[match(cgt.gmf.common,rownames(tdi.emp.df)),]
cgt.sds.df2 = Y.sds.vec[match(cgt.gmf.common,rownames(Y.sds.vec)),]
cgt.tmtA.df2 = Y.tmtA.vec[match(cgt.gmf.common,rownames(Y.tmtA.vec)),]
cgt.tmtB.df2 = Y.tmtB.vec[match(cgt.gmf.common,rownames(Y.tmtB.vec)),]
cgt.fit.df2 = Y.fit.vec[match(cgt.gmf.common,rownames(Y.fit.vec)),]
cgt.dst.df2 = Y.dst.vec[match(cgt.gmf.common,rownames(Y.dst.vec)),]
cgt.pmt.df2 = Y.pmt.vec[match(cgt.gmf.common,rownames(Y.pmt.vec)),]


cgt.gmf.df = data.frame(cbind(Y.tdi.vec.df2,
                              tdi.ncars.df2$ncars,
                              tdi.nhh.df2$nhh,
                              tdi.rent.df2$rent,
                              tdi.emp.df2$emp,
                              cgt.sds.df2,
                              cgt.tmtA.df2,
                              cgt.tmtB.df2,
                              cgt.fit.df2,
                              cgt.dst.df2,
                              cgt.pmt.df2))
colnames(cgt.gmf.df) = c("tdi","ncars","nhh","rent","emp","sds","tmtA","tmtB","fit","dst","pmt")

tdi.cgt.df = data.frame(cbind(Y.tdi.vec.df2,
                              cgt.sds.df2,
                              cgt.tmtA.df2,
                              cgt.tmtB.df2,
                              cgt.fit.df2,
                              cgt.dst.df2,
                              cgt.pmt.df2))
colnames(cgt.gmf.df) = c("tdi","sds","tmtA","tmtB","fit","dst","pmt")

tdi.cgt.df.norm = apply(tdi.cgt.df,2,function(y){
                            norm.vec = y/sqrt(sum(y^2,na.rm=T))
                            return(norm.vec)
                    })

tdi.cgt.df.norm.cov = cov(tdi.cgt.df.norm,use='complete.obs')

tdi.cgt.melt = melt(tdi.cgt.df.norm.cov)
ggplot(data=tdi.cgt.melt,aes(Var2,Var1,fill=value))+
    geom_tile(color="white")+
    scale_fill_gradient(low="blue",high="red",
                        mid="white",midpoint = 0, 
                        limit = c(-1,1), 
                        space = "Lab",
                        name="Pearson\nCorrelation"
                        )+theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1,size = 12, hjust = 1))+coord_fixed()

