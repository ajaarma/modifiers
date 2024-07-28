loc4_0="/hps/software/users/dunham/R_lib/4.0.3"
require('tzdb',lib=loc4_0)
require('readr',lib=loc4_0)
require('ukbtools',lib=loc4_0)
require('data.table',lib=loc4_0)

getUKBB_ICD10 <- function(x) {
    ukb_id = x[1]
    s_name = strsplit(strsplit(x[2],"/")[[1]][11],"_")[[1]][1]
    
    tryCatch(
        {
            df_icd = ukb_icd_diagnosis(my_ukb_data, id = s_name, icd.version = 10)
            s_name_icd = paste(df_icd[,3],collapse="|")
            res = list(ukb_id,s_name,s_name_icd)
            return(res)
        },
        error=function(e){
            message('No ICD10 for ukb: ',ukb_id)
            s_name_icd = "NA"
            res = list(ukb_id,s_name,s_name_icd)
            #print(s_name_icd)
            return(res)
        }
    )
}

# Loading UKBB data
ukb_path="/nfs/research/dunham/samples/ukbb/soft/"
#my_ukb_data = ukb_df("ukb50042.tab",path=ukb_path)

batch_str = gsub("^","batch",seq(1:10)[2:10])

for(batchNum in batch_str) {
    message('Processing Batch Num: ',batchNum)

    outFile = paste0("/nfs/research/dunham/samples/ukbb/data/icd10/",batchNum,"_map.txt",sep="")
    batch1 = paste0("/nfs/research/dunham/samples/ukbb/data/batches/",batchNum,sep="")
    batchMap = "/nfs/research/dunham/samples/ukbb/data/23161_ren_mapFile.txt"

    # Extracting batches and sample names
    bdf_1 = fread(batch1,sep="\t",header=F,stringsAsFactors=F,quote="")
    bmdf = fread(batchMap,sep="\t",header=F,stringsAsFactors=F,quote="")
    bmdf_1 = subset(bmdf,bmdf$V1 %in% as.vector(bdf_1$V1))

    # Loading UKBB data
    ukb_path="/nfs/research/dunham/samples/ukbb/soft/"
    #my_ukb_data = ukb_df("ukb50042.tab",path=ukb_path)

    # Extract ICD10 term per samples
    sample_icd10 = apply(data.frame(bmdf_1),
                         1,function(x) {
                                res = getUKBB_ICD10(x)
                                return(res)
                                }
                        )

    #res.df = data.frame(sample_icd10[[1]],sample_icd10[[2]],sample_icd10[[3]])
    message("Writing out mapped file for: ",batchNum)
    for (i in 1:length(sample_icd10)){
        message("Writing out mapped file for: ",batchNum)
        cat(as.character(sample_icd10[[i]][[1]]),'\t',sample_icd10[[i]][[2]],'\t',sample_icd10[[i]][[3]],'\n',
            file=outFile,append=T)
    }
}
