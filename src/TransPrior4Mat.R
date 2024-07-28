################################################################
# Script to generate matrix for filtered variants Mild &severe
#
# Author: Ajay A. Kumar (EMBL-EBI)
#         aak@ebi.ac.uk
#
###############################################################

args=(commandArgs(TRUE))
loc_4 = "/hps/software/users/dunham/R_lib/4.0.3/"

require('data.table',lib.loc=loc_4)
require('reshape2',lib.loc=loc_4)

`%ni%` <- Negate(`%in%`)

checkIn <- function(x,g,s){
    y = c()
    for(i in x){
        i = as.character(i)
        y = append(y,all(grepl(g,paste0(strsplit(i,s)[[1]]),perl=TRUE)))
    }
    return(y)
}

anal_type = args[[1]]
chr_num = args[[2]]

inp_dir = paste0('/hps/nobackup/dunham/ddd/filterSNV/',anal_type,'/id_genes_burden_V4/',
                 chr_num,'/',collapse="")
het_file_list = list.files(inp_dir,pattern='het.txt')

gene.df = data.frame()
trans_map = '/nfs/research/dunham/resources/ensembl/grch38/ensBioMart_grch38_v100_ENST_lengths_200510.txt'
transData = fread(trans_map,sep='\t',header=T,stringsAsFactors=F,quote="")
transData = data.frame(transData)
colnames(transData) = c("Gene_StableID","Transcript_StableID","Transcript_length",
                        "Transcript_type","Gene_name","Gene_type")

out_file = paste0('/nfs/research/dunham/samples/ddd/analysis/ml_mat/',anal_type,'/',
                 chr_num,'_',anal_type,'_matOutput.csv',collapse="")

log_file = paste0('/nfs/research/dunham/samples/ddd/analysis/ml_mat/',anal_type,'/.',
                 chr_num,'_',anal_type,'_logFile.txt',collapse="")

sink(log_file)

for(ele in het_file_list){ 

    print(ele)
    gene_name = strsplit(ele,'[_]')[[1]][1]
    inp_vars = paste(inp_dir,ele,sep='/') #args[[1]]

    varData = fread(inp_vars,sep='\t',header=T,stringsAsFactors=F,quote="")
    varData = data.frame(varData)

    for (i in 1:dim(varData)[1]){

        trans_ensg_list = strsplit(varData[i,'ANNO_Feature'],'[|]')[[1]]
        trans_pro_cod_index = which(grepl('protein_coding',strsplit(varData[i,'ANNO_BIOTYPE'],'[|]')[[1]])==TRUE)
        trans_impact_index = which(grepl('^HIGH|^MODERATE|^LOW',strsplit(varData[i,'ANNO_IMPACT'],'[|]')[[1]])==TRUE)
        #trans_cano_index = which(grepl('^YES',strsplit(varData[i,'ANNO_CANONICAL'],'[|]')[[1]])==TRUE)
        #pcod_imp_index = intersect(trans_pro_cod_index,intersect(trans_impact_index,trans_cano_index))
        pcod_imp_index = intersect(trans_pro_cod_index,trans_impact_index)

        if (length(pcod_imp_index)!=0){
            trans_pcod_imp_ens_list = trans_ensg_list[pcod_imp_index]
            trans_prior_len_df = transData[transData$Transcript_StableID %in% trans_pcod_imp_ens_list,]
            if (dim(trans_prior_len_df)[1]!=0) {
                trans_prior_sort = trans_prior_len_df[order(trans_prior_len_df$Transcript_length),]
                trans_prior_ens_id = trans_prior_sort$Transcript_StableID[1]
                trans_prior_index = which(trans_ensg_list==trans_prior_ens_id)

                #Extract Transcript related attributes
                gene_transcript = trans_ensg_list[trans_prior_index]
                gene_anno_cons = strsplit(varData[i,'ANNO_Consequence'],'[|]')[[1]][trans_prior_index]
                gene_anno_cadd = strsplit(varData[i,'ANNO_CADD_PHRED'],'[|]')[[1]][trans_prior_index]
                gene_anno_exacPLI = strsplit(varData[i,'ANNO_ExACpLI'],'[|]')[[1]][trans_prior_index]
                gene_anno_revel = strsplit(varData[i,'ANNO_REVEL'],'[|]')[[1]][trans_prior_index]

                #Extract variant related attributes
                var_gnomad = varData[i,'GNOMADgV3_AF']
                var_gnomad_g =  varData[i,'GNOMADg_AF']
                var_gnomad_e =  varData[i,'GNOMADe_AF']
                var_exac =  varData[i,'EXAC_AF']

                #Sample Names
                nhet_list = strsplit(varData[i,'sHet_names_RDAB'],'[;]')[[1]]
                
                #Attribute names
                attr_names = c("Transcript","Consequence","CADD","exacPLI","revel",
                              "gnomADV3_af","gnomAD_g_af","gnomAD_e_af","exac_af"
                              )

                var_id = paste(gene_name,varData[i,'CHROM'],varData[i,'POS'],varData[i,'REF'],
                                varData[i,'ALT'],sep="_")

                attr_names_new = gsub("^",paste0(var_id,'_',collapse=''),attr_names)
                attr_name_val = c(gene_transcript,gene_anno_cons,gene_anno_cadd,
                                  gene_anno_exacPLI,gene_anno_revel,var_gnomad,
                                  var_gnomad_g,var_gnomad_e,var_exac)

                tmp.df = data.frame()
                for (s_name in nhet_list) {
                    tmp_s_list = rep(s_name,length(attr_name_val))
                    tmp.df = data.frame(tmp_s_list,attr_names_new,attr_name_val)
                }
                gene.df = rbind(gene.df,tmp.df)
            }else{
                cat('No Match Transcript in DB: ',gene_name,' ',trans_pcod_imp_ens_list,'\n')
            }
        }else{
            cat('No Transcript found for gene: ',gene_name,'\n')
        }

    }
        
}

geneMat = acast(gene.df,tmp_s_list~attr_names_new,value.var="attr_name_val") 

write.table(geneMat,sep='\t',file=out_file,row.names=T,col.names=T,quote=F)
sink()
