getTransPrior <- function(tmp.var,transData,gene_name) {

    # Subroutine to get Prioritized Transcripts; 
    # Add sep columns & index of Prioritized transcripts

        #message('Get prioritized transcripts')

        transPriorList = c()
        consqPriorList = c()
        hgvspPriorList = c()
        hgvscPriorList = c()
        revelPriorList = c()
        pliPriorList = c()
        amPriorList = c()
        clnPriorList = c()

        varData = data.frame(tmp.var)
        
        for (i in 1:dim(varData)[1]){
            gene_symb_list = strsplit(as.vector(varData[i,'ANNO_SYMBOL']),'[|]')[[1]]
            trans_ensg_list = strsplit(as.vector(varData[i,'ANNO_Feature']),'[|]')[[1]]

            gene_symb_index = which(grepl(gene_name,strsplit(as.vector(varData[i,'ANNO_SYMBOL']),'[|]')[[1]],fixed=T)==TRUE)
            trans_pro_cod_index = which(grepl('protein_coding',strsplit(as.vector(varData[i,'ANNO_BIOTYPE']),'[|]')[[1]])==TRUE)
            trans_impact_index = which(grepl('^HIGH|^MODERATE|^LOW',strsplit(as.vector(varData[i,'ANNO_IMPACT']),'[|]')[[1]])==TRUE)
            trans_cano_index = which(grepl('^YES',strsplit(varData[i,'ANNO_CANONICAL'],'[|]')[[1]])==TRUE)
            #pcod_imp_index = intersect(trans_pro_cod_index,intersect(trans_impact_index,trans_cano_index))
            #pcod_imp_index = intersect(gene_symb_index,intersect(trans_pro_cod_index,trans_impact_index))
            pcod_imp_index = intersect(gene_symb_index,intersect(trans_pro_cod_index,intersect(trans_impact_index,trans_cano_index)))

            
            if (length(pcod_imp_index)!=0){
                trans_pcod_imp_ens_list = trans_ensg_list[pcod_imp_index]
                trans_prior_len_df = transData[transData$Transcript_StableID %in% trans_pcod_imp_ens_list,]
                if (dim(trans_prior_len_df)[1]!=0) {
                    trans_prior_sort = trans_prior_len_df[order(trans_prior_len_df$Transcript_length,decreasing=T),]
                    trans_prior_ens_id = trans_prior_sort$Transcript_StableID[1]
                    trans_prior_index = which(trans_ensg_list==trans_prior_ens_id)

                    #Extract Transcript related attributes
                    gene_transcript = trans_ensg_list[trans_prior_index]
                    transPriorList = c(transPriorList,gene_transcript)

                    #Extract Transcript consequences 
                    gene_anno_cons = strsplit(as.vector(varData[i,'ANNO_Consequence']),'[|]')[[1]][trans_prior_index]
                    consqPriorList = c(consqPriorList,gene_anno_cons)

                    #Extract HGVSp
                    gene_hgvsp = strsplit(as.vector(varData[i,'ANNO_HGVSp']),'[|]')[[1]][trans_prior_index]
                    hgvspPriorList = c(hgvspPriorList,gene_hgvsp)

                    #Extract HGVSc
                    gene_hgvsc = strsplit(as.vector(varData[i,'ANNO_HGVSc']),'[|]')[[1]][trans_prior_index]
                    hgvscPriorList = c(hgvscPriorList,gene_hgvsc)

                    #Extract REVEL score
                    gene_revel = strsplit(as.vector(varData[i,'ANNO_REVEL']),'[|]')[[1]][trans_prior_index]
                    revelPriorList = c(revelPriorList,gene_revel)

                    #Extract pLI score
                    gene_pli = strsplit(as.vector(varData[i,'ANNO_pLI_gene_value']),'[|]')[[1]][trans_prior_index]
                    pliPriorList = c(pliPriorList,gene_pli)

                    #Extract AM score
                    gene_am = strsplit(as.vector(varData[i,'ANNO_am_class']),'[|]')[[1]][trans_prior_index]
                    amPriorList = c(clnPriorList,gene_am)

                    #Extract CLIN score
                    gene_cln = strsplit(as.vector(varData[i,'ANNO_CLIN_SIG']),'[|]')[[1]][trans_prior_index]
                    clnPriorList = c(clnPriorList,gene_cln)
                }else{
                    cat('No Match Transcript in DB: ',gene_name,' ',trans_pcod_imp_ens_list,'\n')
                    gene_transcript = "NA"
                    gene_anno_cons = "NA"
                    gene_hgvsp = "NA"
                    gene_hgvsc = "NA"
                    gene_revel = "NA"
                    gene_pli = "NA"
                    gene_am = "NA"
                    gene_cln = "NA"
                    transPriorList = c(transPriorList,gene_transcript)
                    consqPriorList = c(consqPriorList,gene_anno_cons)
                    hgvspPriorList = c(hgvspPriorList,gene_hgvsp)
                    hgvscPriorList = c(hgvscPriorList,gene_hgvsc)
                    revelPriorList = c(revelPriorList,gene_revel)
                    pliPriorList = c(pliPriorList,gene_pli)
                    amPriorList = c(amPriorList,gene_am)
                    clnPriorList = c(clnPriorList,gene_cln)
                }
            }else{
                cat('No Transcript found for gene: ',gene_name,'\n')
                gene_anno_cons = "NA"
                gene_transcript = "NA"
                gene_anno_cons = "NA"
                gene_hgvsp = "NA"
                gene_hgvsc = "NA"
                gene_revel = "NA"
                gene_pli = "NA"
                gene_am = "NA"
                gene_cln = "NA"
                transPriorList = c(transPriorList,gene_transcript)
                consqPriorList = c(consqPriorList,gene_anno_cons)
                hgvspPriorList = c(hgvspPriorList,gene_hgvsp)
                hgvscPriorList = c(hgvscPriorList,gene_hgvsc)
                revelPriorList = c(revelPriorList,gene_revel)
                pliPriorList = c(pliPriorList,gene_pli)
                amPriorList = c(amPriorList,gene_am)
                clnPriorList = c(clnPriorList,gene_cln)
             }
        }
        
        format_index = which(colnames(varData)=='FORMAT')
        varData_1 = varData[,c(1:format_index-1)]
        varData_2 = varData[,c(format_index:dim(varData)[2])]
        varData.new = cbind(varData_1,transPriorList,consqPriorList,hgvspPriorList,
                            hgvscPriorList,revelPriorList,pliPriorList,amPriorList,clnPriorList,varData_2)

        #return(gene_anno_cons)
        return(varData.new)
    }
