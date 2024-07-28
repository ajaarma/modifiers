#----------------------------------------------------------------------
#Description: Subroutine to filter out top SNPs pertaining to EA
# Filtering criteria:
# (a) compute correlation with NEDU years residualized
# (a) Betas >=0.02 (+ve contribution for Education Attainment
#----------------------------------------------------------------------

require('dplyr')
require('data.table')

betasFile = '/nfs/research/dunham/samples/ukbb/data/gwas/sum_stat/okbay_2022/Betas.csv'
snpsFile = '/hps/nobackup/dunham/ai-uk-can/obkay_2022/plink_conv/maf_005/snpQC.snplist'
varGeneFile = '/hps/nobackup/dunham/ai-uk-can/obkay_2022/obk_vars_gene_map.txt'
geneVarMap = '/hps/nobackup/dunham/ai-uk-can/obkay_2022/obk_gene_var_map_sort.txt'

betaDat = read.table(betasFile,sep=",",header=T,stringsAsFactors=F)
snpsQC = read.table(snpsFile,sep=":",header=F,stringsAsFactors=F)
snpsGene = read.table(varGeneFile,sep='\t',header=F,stringsAsFactors=F)
geneSnps = read.table(geneVarMap,sep='\t',header=F,stringsAsFactors=F)

betaDat$chr_pos_var = paste(betaDat$chr_name,betaDat$chr_position,betaDat$effect_allele,betaDat$noneffect_allele,sep="_")
betaDat$chr_pos_var = gsub("chr","",betaDat$chr_pos_var)

betaDat$chr_pos = paste(betaDat$chr_name,betaDat$chr_position,sep="_")
betaDat$chr_pos = gsub("chr","",betaDat$chr_pos)

colnames(snpsQC) = c('chr_num','pos_var')
snpsQC$chr_pos_var = paste(snpsQC$chr_num,snpsQC$pos_var,sep="_")
snpsQC$pos = apply(data.frame(snpsQC$pos_var),1,function(x){tmp = strsplit(x,"_")[[1]];return(tmp[1])})
snpsQC$chr_pos = paste(snpsQC$chr_num,snpsQC$pos,sep="_")

colnames(snpsGene) = c('chr_pos_var','gene_names')
colnames(geneSnps) = c('gene_names','vars','count')

#-----------------------------------
# Commonality between above dataset
#-----------------------------------

beta_snpsqc = dplyr::inner_join(betaDat,snpsQC,by=join_by(chr_pos_var == chr_pos_var))
colnames(beta_snpsqc)[8] = 'chr_pos'
snpsqc_diff = snpsQC[ !snpsQC$chr_pos_var %in% betaDat$chr_pos_var,]
beta_snpsqc_diff = dplyr::inner_join(betaDat,snpsqc_diff,by=join_by(chr_pos == chr_pos))

#----------------------------------
# 3811 common variant (MAF >5%)
#-----------------------------------
beta.snpsqc.merge = unique(rbind(beta_snpsqc[,c(1:8)],beta_snpsqc_diff[,c(1:8)]))


#------------------------------------------
# Criteria 1: Select SNPs with large effect
#------------------------------------------
a1.cv = beta.snpsqc.merge[order(beta.snpsqc.merge$eaf,decreasing=T),]

#--------------------------------------------------------------------------
# Criteria 2: Select SNPs having high degree of correlation with NEDU years
#--------------------------------------------------------------------------

#rv_cv_list = '/hps/nobackup/dunham/ai-uk-can/obkay_2022/dl_data/cv_rvframe_nedu_list'
#load(rv_cv_list)

#edu.cv.rv = cv_rv_nedu_list[[1]]
#nedu.rsid = cv_rv_nedu_list[[2]]

#edu.cv.rv.cor = apply(data.frame(edu.cv.rv[,-1]),2,function(x){cor_val = cor(as.numeric(x),nedu.rsid[,3]);return(cor_val)})

colnames(edu.cv.rv)[5644] = 'XNF1_RV'
edu.cv.rv.cor = c()

for(i in 2:dim(edu.cv.rv)[2]){
    print(i)
    var_gt = as.numeric(edu.cv.rv[,i])
    #var_gt[var_gt==2] = 1 #Gives 70 snps in total
    #cor_val = cor(as.numeric(edu.cv.rv[,i]),nedu.rsid[,3])
    cor_val = cor(as.numeric(var_gt),nedu.rsid[,3])
    edu.cv.rv.cor = c(edu.cv.rv.cor,cor_val)
}

edu.var.cor = data.frame(cbind(colnames(edu.cv.rv)[-1],edu.cv.rv.cor))
colnames(edu.var.cor) = c('vars','cor_value')
edu.var.cor$cor_value = as.numeric(edu.var.cor$cor_value)
edu.var.cor$vars = gsub("X","",edu.var.cor$vars)

edu.var.cor$chr_pos = apply(data.frame(edu.var.cor$vars),1,function(x){res = strsplit(x,"_")[[1]];out = paste(res[1],res[2],sep="_");return(out)})
edu.var.cor.uniq = unique(edu.var.cor)
edu.var.cor.uniq.sort = edu.var.cor.uniq[order(edu.var.cor.uniq$cor_value,decreasing=T),]

a1.cv.cor = dplyr::inner_join(a1.cv,edu.var.cor.uniq.sort,by=join_by(chr_pos == chr_pos))


edu.top.snp = subset(a1.cv.cor,cor_value >=0.01 | cor_value <=-0.01)

#-------------------------------------------
# Map top correlated snps to Ensembl GeneIDs
#-------------------------------------------
snpsGene$chr_pos = apply(data.frame(snpsGene$chr_pos_var),1,function(x){res = strsplit(x,"_")[[1]];out = paste(res[1],res[2],sep="_");return(out)})
edu.top.snp.gene = dplyr::inner_join(edu.top.snp,snpsGene,by=join_by(chr_pos == chr_pos))

#--------------------------------------------
# Get mapped Genes and snps
#-------------------------------------------

nedu.top.snp.geneList = c()
for (i in 1:dim(edu.top.snp)[1]){
    message('Mappedd gene num: ',i)
    strs = strsplit(edu.top.snp.gene$gene_names[i],",")[[1]]
    for (e in strs){
        if(e !=''){
            e = gsub(" ","",e)
            nedu.top.snp.geneList = c(nedu.top.snp.geneList,e)
        }
    }
}

top.snp.gene.df = data.frame(nedu.top.snp.geneList)
colnames(top.snp.gene.df) = "gene_names"

top.snp.gene.vars = dplyr::left_join(top.snp.gene.df,geneSnps,by=join_by(gene_names == gene_names))

out_nedu_dir = '/hps/nobackup/dunham/ai-uk-can/obkay_2022/dl_data/nedu'
top.gene.file = paste(out_nedu_dir,'top.snp.gene.map.txt',sep="/")

top.snp.gene.vars.sort = top.snp.gene.vars[order(top.snp.gene.vars$count,decreasing=T),]
write.table(top.snp.gene.vars.sort,file=top.gene.file,row.names=F,col.names=T,quote=F)

# Extract top SNPs and Genes and create genotype matrix
top.gene.gt.mat = getGenotype(edu.cv.rv,top.snp.gene.vars)
write.table(top.gene.gt.mat,sep=',',file=paste(out_nedu_dir,'genotype.csv',sep='/'),row.names=F,col.names=F,quote=F)

# Write the snps size file
geneNumList = top.snp.gene.vars$count
geneNumList = c(geneNumList,1)
write.table(geneNumList,file=paste(out_nedu_dir,'snp_size.csv',sep='/'),row.names=F,col.names=F,quote=F)
#Write all filtered top correlated SNPs data with mappings
write.table(edu.top.snp.gene,sep='\t',file=paste(out_nedu_dir,'edu.top.snp.gene.txt',sep='/'),row.names=F,col.names=T,quote=F)


neduList = list()
neduList[[1]] = top.snp.gene.vars.sort
neduList[[2]] = top.gene.gt.mat
neduList[[3]] = edu.top.snp.gene

save(neduList,file=paste(out_nedu_dir,'neduList',sep='/'))

getGenotype <- function(edu.cv.rv,top.snp.gene.vars){

    gt.df = data.frame(edu.cv.rv[,1])
    gt_col_names = c("f.eid_samples")

    for (i in 1:dim(top.snp.gene.vars)[1]){
        vars_strs = strsplit(top.snp.gene.vars$vars[i],",")[[1]]
        for(vars in vars_strs){
            ind = which(grepl(vars,colnames(edu.cv.rv))==TRUE)[1]
            gt.df = cbind(gt.df,edu.cv.rv[,ind])
            gt_col_names = c(gt_col_names,vars)
        }
    }
    nf1.ind = 5644
    gt.df = cbind(gt.df,edu.cv.rv[,nf1.ind])
    gt_col_names = c(gt_col_names,colnames(edu.cv.rv)[5644])

    colnames(gt.df) = gt_col_names

    return(gt.df[,-1])
}
