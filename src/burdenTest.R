lib_loc ="/hps/software/users/dunham/R_lib/4.0.3/" 
#require('data.table',lib.loc="/hps/software/users/dunham/R_lib/4.0.3/")
require("fuzzySim",lib=lib_loc)

#geneList = "/nfs/research/dunham/samples/ddd/data/gene_list/id_high_cap.txt"
#id_hc = read.table(geneList,sep="\t",header=T)

#mild_severe_list = c("gdd_mild","gdd_severe")

#for (ms_type in mild_severe_lis) {
#    ms_path = paste("/hps/nobackup/dunham/ddd/filterSNV",ms_type,"merge_data/",sep="/")
#    ms_abs = strsplit(ms_type,"_")[[1]][2]
#    ms_file = paste(ms_path,
#                    paste("merged.all.AC0.norm.anno.sites.gdd.freq.",
#                    ms_abs,".exon.0_01.imp.prior.tab",sep=""),
#                    sep="")
#    raw_variants <- fread(ms_file, header = T, stringsAsFactors = F, quote="")

#}

msFile = "/hps/nobackup/dunham/ddd/filterSNV/gene_count_full_list_burden/ms_gene_var_count.txt"
msCount = read.table(msFile,sep="\t")
colnames(msCount) = c("ChrNum","Gene","Mild","Severe")

resChi = c()
resFe = c()
for (i in 1:dim(msCount)[1]){
    
    mat = mat.or.vec(2,2)
    rownames(mat) = c("gdd_mild","gdd_severe")
    colnames(mat) = c("CountM","TotalSamples")
    mat[1,1] = msCount[i,3]
    mat[1,2] = 639
    mat[2,1] = msCount[i,4]
    mat[2,2] = 1038

    print(mat)
    res.chi = chisq.test(mat,correct=F,simulate.p.value=TRUE)
    res.fe = fisher.test(mat,simulate.p.value=T)
    p_val_chi = res.chi$p.value
    p_val_fe = res.fe$p.value
    resChi = c(resChi,p_val_chi)
    resFe = c(resFe,p_val_fe)
}

msResult = cbind(msCount,resChi,resFe)
p.adj.val = p.adjust(as.vector(msResult[,6]),method=)
msResult = cbind(msResult,p.adj.val)
colnames(msResult)[5] = "p_Value_chi"
colnames(msResult)[6] = "p_Value_fe"
colnames(msResult)[7] = "p_value_fe_adj"

outFile = "/hps/nobackup/dunham/ddd/filterSNV/gene_count_full_list_burden/ms_gene_var_count_sig.txt"
write.table(msResult,sep="\t",file=outFile,row.names=F,col.names=T,quote=F)

out_fdr = "/hps/nobackup/dunham/ddd/filterSNV/gene_count_full_list_burden/ms_gene_var_count_sig_fdr.txt"
df_05 = msResult[msResult$p_Value_fe < 0.051,]
df_05_pval = df_05$p_Value_fe
fdr_pval = p.adjust(df_05_pval,method="fdr")

df_05_fdr = cbind(df_05,fdr_pval)
write.table(df_05_fdr,sep="\t",file = out_fdr,row.names=F,col.names=T,quote=F)

