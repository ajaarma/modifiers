###################################################################
#
# Description: Script to process ID & IQ gene list. 
# Get relevant unique,shared list of genes for both ID and IQ
#
# Author: aak@ebi.ac.uk
#
#
###################################################################

require('ggvenn')
require('ggplot2')


id_genes = read.table('/nfs/research/dunham/samples/ddd/data/gene_list/id_mr_list_chrMap_prune.txt',sep='\t')
iq_genes = read.table('/nfs/research/dunham/samples/ddd/data/gene_list/intel_gene_list_chrMapSort.txt',sep='\t')

id_uniq = setdiff(id_genes$V1, iq_genes$V1)
iq_uniq = setdiff(iq_genes$V1, id_genes$V1)
id_iq_genes = intersect(id_genes$V1,iq_genes$V1)

id_genes_uniq = id_genes[id_genes$V1 %in% id_uniq,]
iq_genes_uniq = iq_genes[iq_genes$V1 %in% iq_uniq,]
id_iq_genes_shared = id_genes[id_genes$V1 %in% id_iq_genes,]


write.table(id_genes_uniq,sep='\t',file='/nfs/research/dunham/samples/ddd/data/gene_list/ID_IQ/id_genes_uniq.txt',col.names=F,row.names=F,quote=F)
write.table(iq_genes_uniq,sep='\t',file='/nfs/research/dunham/samples/ddd/data/gene_list/ID_IQ/iq_genes_uniq.txt',col.names=F,row.names=F,quote=F)
write.table(id_iq_genes_shared,sep='\t',file='/nfs/research/dunham/samples/ddd/data/gene_list/ID_IQ/id_iq_genes_shared.txt',col.names=F,row.names=F,quote=F)

list_venn = list(ID=sort(id_genes$V1),IQ=sort(iq_genes$V1))
ggsave('/nfs/research/dunham/samples/ddd/data/gene_list/ID_IQ/id_iq_genes.png')
ggvenn(list_venn,c("ID","IQ"))
