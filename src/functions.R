require('data.table')
require('dplyr')

getClinvar1 <- function(df.gene){
        
    cat('--Inside Clinvar 1 --','\n')
    if (!"ANNO_CLIN_SIG" %in% colnames(df.gene)) {
        cat("  --ANNO_CLIN_SIG anno not found --",'\n')

        df.gene.clinvar.path.1 = data.frame()
        return(df.gene.clinvar.path.1)
    }
  
	#Extract Pathogenic/Likely Patho from Anno-Clinsig
	df.gene.clinvar.path.1 = df.gene[apply(data.frame(df.gene[,c("ANNO_SYMBOL","ANNO_CLIN_SIG")]),1,function(x){
                                                 gene_list = x[1]
                                                 clin_sig_list = x[2]
                                                 gene_index = grep(gene_id.app,strsplit(gene_list,"[|]")[[1]])
                                                 clin_sig_strs = strsplit(clin_sig_list,"[|]")[[1]][gene_index]
                                                 clin_sig_strs_and = unlist(apply(data.frame(clin_sig_strs),
                                                                                  1,function(x){
                                                                                      y = strsplit(x,"[&]")[[1]]
                                                                                    }
                                                                                  )
                                                                           )
                                                 y = any(grepl('pathogenic$|likely_pathogenic$',clin_sig_strs_and))
                                                 return(y)
                                                 }
                                             ),]

	return(df.gene.clinvar.path.1)	

}

getClinvar2 <- function(df.gene){
        
    cat('--Inside Clinvar 2 --','\n')
    if (!"CLNSIG" %in% colnames(df.gene)) {
        cat("  --CLNSIG anno not found --",'\n')

        df.gene.clinvar.path.2 = data.frame()
        return(df.gene.clinvar.path.2)
    }
    
	df.gene.clinvar.path.2 = df.gene[apply(data.frame(df.gene[,c("ANNO_SYMBOL","CLNSIG")]),1,function(x){
										gene_list = x[1]
										clin_sig_list = x[2]
										gene_index = grep(gene_id.app,strsplit(gene_list,"[|]")[[1]])
										clin_sig_strs = strsplit(clin_sig_list,"[|]")[[1]][gene_index]
										clin_sig_strs_and = unlist(apply(data.frame(clin_sig_strs),
																			1,function(x){
																			y = strsplit(x,"[/]")[[1]]
																		}
																	)
																)
										y = any(grepl('Pathogenic|Likely_pathogenic',clin_sig_strs_and))
										return(y)
										}
								),]

	return(df.gene.clinvar.path.2)
	}

getClinvar3 <- function(df.gene) {

        cat('--Inside Clinvar 3 --','\n')
        if (!"CLNSIGCONF" %in% colnames(df.gene)) {
            cat("  --CLNSIGCONF anno not found --",'\n')
            df.gene.clinvar.path.3 = data.frame()
            return(df.gene.clinvar.path.3)
        }
		df.gene.clinvar.path.3 = df.gene[apply(data.frame(df.gene[,c("ANNO_SYMBOL","CLNSIGCONF")]),1,function(x){
											gene_list = x[1]
											clin_sig_list = x[2]
											gene_index = grep(gene_id.app,strsplit(gene_list,"[|]")[[1]])
											clin_sig_strs = strsplit(clin_sig_list,"[|]")[[1]][gene_index]
											clin_sig_strs_and = unlist(apply(data.frame(clin_sig_strs),
																				1,function(x){
																					y = strsplit(x,"[,]")[[1]]
																					}
																			)
																	)
											y = any(grepl('Pathogenic|Likely_pathogenic',clin_sig_strs_and))
											return(y)
												}
										),]


		return(df.gene.clinvar.path.3)
	}

#Added functionality of filtering by AlphaM pathogenicity class
#11.11.2023
getAlphaM <- function(df.gene) {

        cat("--Inside AlphaM --",'\n')
        if (!"ANNO_am_class" %in% colnames(df.gene)) {
            cat("  --AlphaM anno not found --",'\n')
            
            df.gene.am.path = data.frame()
            return(df.gene.am.path)
        }
		df.gene.am.path = df.gene[apply(data.frame(df.gene[,c("ANNO_SYMBOL","ANNO_am_class")]),1,function(x){
											gene_list = x[1]
											am_sig_list = x[2]
											gene_index = grep(gene_id.app,strsplit(gene_list,"[|]")[[1]])
											am_sig_strs = strsplit(am_sig_list,"[|]")[[1]][gene_index]
											am_sig_strs_and = unlist(apply(data.frame(am_sig_strs),
																				1,function(x){
																					y = strsplit(x,"[,]")[[1]]
																					}
																			)
																	)
											y = any(grepl('pathogenic|likely_pathogenic',am_sig_strs_and))
											return(y)
												}
										),]

		return(df.gene.am.path)
	}


getCadd <- function(df.gene.path){

	df.gene.cadd = df.gene.path[apply(data.frame(df.gene.path[,c("ANNO_SYMBOL","ANNO_CADD_PHRED")]),
									1,function(x){
										gene_list = x[1]
										cadd_list = x[2]
										gene_index = grep(gene_id.app,strsplit(gene_list,"[|]")[[1]])
										cadd_gene_strs = strsplit(cadd_list,"[|]")[[1]][gene_index]
										y = any(as.numeric(cadd_gene_strs) >=25.0)
										return(y)
									}
								),]


}

filterRD_AB <- function(sample_gt,minReads,abCutoff){

    # Subroutine to filter by Allelic Balance > 0.15

    cat('--Inside filterRD_AB --','\n')

    ref_ab_list = c()
    alt_ab_list = c()
    gt_list = c()
    s_list = c()

    if (dim(sample_gt)[2]!=0) {
        for(i in 1:dim(sample_gt)[2]){
            gt = sample_gt[1,i]
            s_name = colnames(sample_gt)[i]
            gt_strs = strsplit(as.vector(gt),":")
            ref_alt_strs = strsplit(gt_strs[[1]][3],",")[[1]]
            ref_rd = as.numeric(ref_alt_strs[1])
            alt_rd = as.numeric(ref_alt_strs[2])
            total_rd = as.numeric(gt_strs[[1]][2])
            ref_ab = ref_rd/total_rd
            alt_ab = alt_rd/total_rd

            if (alt_rd >=minReads || alt_ab >=abCutoff){
                gt_list = c(gt_list,gt)
                ref_ab_list = c(ref_ab_list,ref_ab)
                alt_ab_list = c(alt_ab_list,alt_ab)
                s_list = c(s_list,s_name)
            }

        }
    }
    if (length(s_list)==0){
        s_list = c("NA");gt_list=c("NA");ref_ab_list = c("NA");alt_ab_list = c("NA")
    }
    return(c(paste0(s_list,collapse=";"),
             paste0(gt_list,collapse=";"),
             paste0(ref_ab_list,collapse=";"),
             paste0(alt_ab_list,collapse=";")
           )
    )
}

