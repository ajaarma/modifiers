require('ukbtools')
#tab_path = '/nfs/research/dunham/samples/ukbb/soft/basket_2017324'
##path_to_example_data <- system.file("extdata", package = "ukbtools")

##df <- ukb_df("ukb50042", path = path_to_example_data)


######## Get all ICD10 codes and meaning from ukbtools ######
ukb_icd10_file = '/nfs/research/dunham/samples/ukbb/soft/basket_2017324/ukbb_icd10_codes.txt'
icd10.codes = ukbtools::icd10codes
write.table(icd10.codes,sep='\t',file=ukb_icd10_file,row.names=F,col.names=F,quote=F)

### ICD10 terms for all the samples ####
##samples_icd10 = ukb_icd_diagnosis(df1,id = df1$eid, icd.version=10)

samplesUniq = unique(samples_icd10$sample)

#samplePheno.df = data.frame()
#sampleList = c()
#phenoList = c()

#count = 1
#for (s in samplesUniq){
#    cat(s,'\t',count,'\n')
#    pheno = paste(samples_icd10$meaning,collapse="|")
#    sampleList = c(sampleList,s)
#    phenoList = c(phenoList,pheno)
#    count = count+1
#}

#samplePheno.df = cbind(sampleList,phenoList)
#colnames(samplePheno.df) = c("sample","pheno")

#write.table(samplePheno.df,sep='\t',row.names=F,col.names=F,quote=F)

phenoList = apply(data.frame(samplesUniq),1,function(x){
                                pheno = paste(samples_icd10[samples_icd10$sample==x,'meaning'],collapse="|")
                                return(pheno)
                                })

samplePheno.df = cbind(samplesUniq,phenoList)
s_pheno_file = '/nfs/research/dunham/samples/ukbb/soft/basket_2017324/samplePheno.txt'
write.table(samplePheno.df,file=s_pheno_file,sep='\t',row.names=F,col.names=F,quote=F)


