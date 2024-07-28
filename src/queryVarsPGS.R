require('data.table')

#------------------------------------------------------------------
# Description: Query genes and variants in GWAS common variants
#           against GWAS
#
#------------------------------------------------------------------
pgsFile = '/hps/nobackup/dunham/ai-uk-can/obkay_2022/ptv/py_er/GCST90105038.tsv'
ptvFile = '/hps/nobackup/dunham/ai-uk-can/obkay_2022/ptv/py_er/pgs_ptv_genes_py.txt'
gwasData = fread(pgsFile,sep='\t',header=T,stringsAsFactors=T)
ptvData = read.table(ptvFile,header=T,sep='\t')

fitList = list()
for(i in 1:dim(ptvData)[1]){
    geneList = ptvData[i,15]
    geneStrs = gsub(" ","",strsplit(geneList,',')[[1]])
    tmpList = list()
    count = 1
    for (ele in geneStrs){
        tmp = subset(gwasData,grepl(ele,MAPPED_GENE))
        tmpList[[count]] = tmp
        count = count+1
        
    }
    tmpDF = rbindlist(tmpList)
    fitList[[i]] = tmpDF

}

#---------------------------------------------
# Pathway: MAPK signaling pathway
# FIT gene: MAPK1, chr22
#--------------------------------------------

#MAPK signaling pathway
gwas_genes_mapk_py = fitList[[1]]
gwas_gene_in_chr22 = subset(gwas_genes_mapk_py,CHR_ID==22)$MAPPED_GENE #Yields: PDGFB - RPL3
write.table(gwas_genes_mapk_py,sep='\t',file='FIT_EA_GWAS_shared_mapk_pathway.txt')

#Proteoglycans in Cancer
gwas_genes_pg_py = fitList[[4]]
gwas_gene_in_chr22 = subset(gwas_genes_pg_py,CHR_ID==22)$MAPPED_GENE #Yields: SHANK3, WNT7B-LINC00899
write.table(gwas_genes_pg_py,sep='\t',file='FIT_EA_GWAS_shared_pg_pathway.txt')








