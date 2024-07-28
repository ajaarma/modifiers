###############################################################################
# Script to Prioritize Transcripts for the filtered set of variants in UKBB
# Also, generates the bar plots for each of the variant consequences 
# Applicable for both ID and IQ genes
#
# Author: Ajay A. Kumar (EMBL-EBI)
#         aak@ebi.ac.uk
#
###############################################################################

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

getVarDataGene <- function(gene_id,pos,chr_num){

    #message('Inside variant Data per batches')
    for (k in 1:10) {
        anal_type = paste0('ukbb_batch',k,collapse='') #args[[1]]
        run_type = 'id_genes' #args[[2]]

        var_file = paste0('/hps/nobackup/dunham/ukbb/filterSNV/',anal_type,'/',run_type,'/',chr_num,'/',gene_id,'_variants.txt',collapse="")
        if(file.exists(var_file)) {
            #het_file_list = list.files(inp_dir,pattern='_variants.txt')
            varData = fread(var_file,sep='\t',head=T,stringsAsFactors=T)
            varData = data.frame(varData)

            tmp.var = subset(varData,POS==pos)
            if(dim(tmp.var)[1] !=0) {
                trans_consq = getTransPrior(tmp.var,transData,gene_id)
                break
            }else{
                message('Searching in next batch')
                next
            }
        }else {
            next
        }
    }
    return(trans_consq)
}

#write.table(gene.df,sep='\t',file=out_file,row.names=T,col.names=T,quote=F)

getTransPrior <- function(tmp.var,transData,gene_name) {

        #message('Get prioritized transcripts')

        varData = data.frame(tmp.var)
        for (i in 1:dim(varData)[1]){

            trans_ensg_list = strsplit(as.vector(varData[i,'ANNO_Feature']),'[|]')[[1]]
            trans_pro_cod_index = which(grepl('protein_coding',strsplit(as.vector(varData[i,'ANNO_BIOTYPE']),'[|]')[[1]])==TRUE)
            trans_impact_index = which(grepl('^HIGH|^MODERATE|^LOW',strsplit(as.vector(varData[i,'ANNO_IMPACT']),'[|]')[[1]])==TRUE)
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
                    gene_anno_cons = strsplit(as.vector(varData[i,'ANNO_Consequence']),'[|]')[[1]][trans_prior_index]

                }else{
                    cat('No Match Transcript in DB: ',gene_name,' ',trans_pcod_imp_ens_list,'\n')
                    gene_anno_cons = "NA"
                }
            }else{
                cat('No Transcript found for gene: ',gene_name,'\n')
                gene_anno_cons = "NA"
            }
        }
        return(gene_anno_cons)
    }

getColumnPlot <- function(inp.df,case_id,out_plot){
    require('ggplot2')
    # Subroutine to generate the distribution of Consequences 
    # of the variants

    #a11 = table(inp.df$Consequences)
    #a11.df = data.frame(names(a11),as.vector(a11))
    a11.df = inp.df

    #colnames(a11.df) = c('Consq','Values')
    g = ggplot(a11.df,aes(Consq,Prop,fill=Cases))+geom_col(position='dodge2')+
        labs(x="\n Variant Consequences\n",
             y="\n Proportion(%)",
             #title=paste0("Distribution of Variant Consequences for: ",case_id,collapse=""))+
             title=paste0(case_id,collapse=""))+
        theme(plot.title = element_text(face="bold",hjust=0.5,size=18))+
        theme(axis.text.x=element_text(face="bold",angle = -90, hjust = 0,size=11))+
        theme(axis.text.y=element_text(face='bold',angle = 0, size=11))+
        theme(axis.title=element_text(face='bold',size=16))+
        theme(legend.title=element_text(face='bold',size=16))+
        theme(legend.text=element_text(face='bold',size=12))
    ggsave(g,file=out_plot,width=20,height=12,units='in')
    #print(g)
}

run_type = 'id_genes' #args[[1]]
caseList = c('Case1','Case2','Case3')

outCaseDF = data.frame()
for (case_id in caseList) {
    print(case_id)
    #Case IDs
    if(case_id =='Case1'){
        case_id_ext = '_vt_rc25_hc6'
    }else if(case_id =='Case2'){
        case_id_ext = '_vt_urcgte25_hc6'
    }else if(case_id =='Case3'){
        case_id_ext = '_vt_urclte25_hc6'
    }

    #ID/IQ genes input/output directories
    if (run_type =='id_genes'){
        out_dir_ext = 'id_genes_ukbb_merge_all'
        plot_str = 'ID genes, All 3 Cases'
    }else if (run_type =='int_genes'){
        out_dir_ext = 'int_genes_ukbb_merge_all'
        plot_str = 'IQ genes, ALL 3 Cases'
    }

    #Reading Transcript data
    trans_map = '/nfs/research/dunham/resources/ensembl/grch38/ensBioMart_grch38_v100_ENST_lengths_200510.txt'
    transData = fread(trans_map,sep='\t',header=T,stringsAsFactors=F,quote="")
    transData = data.frame(transData)
    colnames(transData) = c("Gene_StableID","Transcript_StableID","Transcript_length",
                        "Transcript_type","Gene_name","Gene_type")

    out_file = paste0('/hps/nobackup/dunham/ukbb/filterSNV/',out_dir_ext,'/',case_id,'_map_variants_consq.txt',collapse="")
    #out_plot = paste0('/hps/nobackup/dunham/ukbb/filterSNV/',out_dir_ext,'/',case_id,'_variants_consq_dist.png',collapse="")
    out_plot = paste0('/hps/nobackup/dunham/ukbb/filterSNV/',out_dir_ext,'/','All_cases_variants_consq_dist.png',collapse="")
    #out_plot = paste0('/hps/nobackup/dunham/ukbb/filterSNV/',out_dir_ext,'/',collapse="")
    #log_file = paste0('/hps/nobackup/dunham/ukbb/filterSNV/tryOuts/',chr_num,'_',anal_type,'_logFile.txt',collapse="")

    #Read the variants
    inpVar = paste0('/hps/nobackup/dunham/ukbb/filterSNV/',out_dir_ext,'/merge_all',case_id_ext,'.txt',collapse='')
    varData = fread(inpVar,sep='\t',header=F,stringsAsFactors=T)
    varList = unique(varData$V1)

    transConsqList = c()
    count = 1
    for (var_id in varList){
        chr_num = strsplit(as.vector(var_id),'[_]')[[1]][1]
        pos = strsplit(as.vector(var_id),'[_]')[[1]][2]
        gene_id = strsplit(as.vector(var_id),'[_]')[[1]][3]

        trans_consq = getVarDataGene(gene_id,pos,chr_num)
        transConsqList = c(transConsqList,trans_consq)
        count = count+1
    }

    caseVec = rep(case_id,length(transConsqList))
    out.df = data.frame(cbind(as.vector(varList),transConsqList))
    colnames(out.df) = c('Variants','Consequences')

    a11 = table(out.df$Consequences)
    a11.frac = 100*a11/length(varList)
    caseVec = rep(case_id,length(a11))
    a11.df = data.frame(names(a11),as.vector(a11),as.vector(a11.frac),caseVec)
    colnames(a11.df) = c('Consq','Freq','Prop','Cases')
    
    outCaseDF = rbind(outCaseDF,a11.df)
    write.table(out.df,file=out_file,sep='\t',col.names=T,row.names=F,quote=F)
    #getColumnPlot(out.df,case_id,out_plot)
}
    getColumnPlot(outCaseDF,plot_str,out_plot)

