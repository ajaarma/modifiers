args=(commandArgs(TRUE))
loc_4 = "/hps/software/users/dunham/R_lib/4.0.3/"

require('data.table',lib.loc=loc_4)

`%ni%` <- Negate(`%in%`)

checkIn <- function(x,g,s){
    y = c()
    for(i in x){
        i = as.character(i)
        y = append(y,all(grepl(g,paste0(strsplit(i,s)[[1]]),perl=TRUE)))
    }
    return(y)
}


#tmp_data = "/hps/nobackup/dunham/ukbb/filterSNV/ukbb_batch1/tmp_data/"
raw_var = args[[1]]
chrNum = args[[2]]
outDir = args[[3]]

#outDir = "/hps/nobackup/dunham/ukbb/filterSNV/ukbb_batch1/id_genes/"
#chrNum=20
#raw_var = paste('/hps/nobackup/dunham/ukbb/filterSNV/ukbb_batch1/tmp_data/',chrNum,
#raw_var = paste(tmp_data,chrNum,
#                '/tmpFile_',chrNum,'.AC0.norm.vep.anno.sites.ukbb.freq.batch1.exon.0_01.imp.tab',sep="")

geneList = '/nfs/research/dunham/samples/ddd/data/gene_list/id_mr_list_chrMap_prune.txt'

varData = fread(raw_var,sep="\t",header=T,stringsAsFactors=F,quote="")
message("--Data Loaded\n")
geneChrom = fread(geneList,sep="\t",header=F,stringsAsFactors=F,quote="")

geneList = subset(geneChrom,grepl(chrNum,geneChrom$V2,fixed=TRUE))
nHetsList = c()
nHomsList = c()

for (i in 1:dim(geneList)[1]){
    gene_name = geneList[i,1]
    message('--GeneName:',gene_name,'\n')
    gene_id = paste("",geneList[i,1],"$",sep="")
    df.gene = subset(varData,grepl(gene_id, varData$ANNO_SYMBOL))
    df.gene.path = subset(df.gene,(grepl("Pathogenic$",df.gene$ANNO_CLIN_SIG,perl=T,ignore.case=T) | 
                                   grepl("Pathogenic$",df.gene$CLNSIG,perl=T,ignore.case=F) |
                                   grepl("Pathogenic",df.gene$CLNSIGCONF,perl=T,ignore.case=F) 
                                   )
                         )

    nHetsList = c()
    nHomsList = c()
    
    sample_index = which(colnames(df.gene.path)=="FORMAT")+1
    df.cols = dim(df.gene.path)[2]
    df.gene.path = as.data.frame(df.gene.path)

    if (dim(df.gene.path)[1]!=0){
        for (j in 1:dim(df.gene.path)[1]){
            
            sampleGT = df.gene.path[j,c(sample_index:df.cols)]
            tmp_nHets = apply(data.frame(sampleGT),1,function(x) grepl('^0/1:|^1/0:|^1:',x,perl=TRUE))
            sHet_names = sampleGT[j,which(tmp_nHets==TRUE)]
            
            nHets = length(which(tmp_nHets==TRUE))
            tmp_nHoms = apply(data.frame(sampleGT),1,function(x) grepl('^0/0:|^1/1:|^0:',x,perl=TRUE))
            nHoms = length(which(tmp_nHoms==TRUE))
            
            nHetsList = c(nHetsList,nHets)
            nHomsList = c(nHomsList,nHoms)
        }

        #df.gene.path["nHets"] = nHetsList
        #df.gene.path["nHoms"] = nHomsList
        df.gene.path = cbind(df.gene.path,nHetsList,nHomsList)
        out_file = paste(outDir,'/',gene_name,'_variants.txt',sep="")
        write.table(df.gene.path,sep="\t",file=out_file,col.names=T,row.names=F,quote=F)
    }
    else{
        message('--No Pathogenic Variants found for:',gene_name,'\n')
    }

}

