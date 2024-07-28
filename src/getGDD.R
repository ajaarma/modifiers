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
#raw_var = "/hps/nobackup/dunham/ddd/filterSNV/gdd_severe/tmp_data/chrX/tmpFile_chrX.AC0.norm.vep.anno.sites.gdd.freq.severe.exon.0_01.imp.tab"
#chrNum = "X"
#outDir = "/hps/nobackup/dunham/ddd/filterSNV/gdd_severe/id_genes/chrX"
raw_var = args[[1]]
chrNum = args[[2]]
outDir = args[[3]]

#outDir = "/hps/nobackup/dunham/ukbb/filterSNV/ukbb_batch1/id_genes/"
#outDir = "/hps/nobackup/dunham/ddd/filterSNV/gdd_severe/id_genes/20"
#chrNum=20
#raw_var = paste('/hps/nobackup/dunham/ukbb/filterSNV/ukbb_batch1/tmp_data/',chrNum,
#raw_var = paste(tmp_data,chrNum,
#                '/tmpFile_',chrNum,'.AC0.norm.vep.anno.sites.ukbb.freq.batch1.exon.0_01.imp.tab',sep="")

filterRD_AB <- function(sample_gt,minReads,abCutoff){
    
    # Subroutine to filter by Allelic Balance > 0.15 
    
    ref_ab_list = c()
    alt_ab_list = c()
    gt_list = c()
    s_list = c()

    if (dim(sample_gt)[2]!=0) {   
        for(i in 1:dim(sample_gt)[2]){
            gt = sample_gt[1,i]
            s_name = colnames(sample_gt)[i]
            gt_strs = strsplit(as.vector(gt),":")
            ref_alt_strs = strsplit(gt_strs[[1]][2],",")[[1]]
            ref_rd = as.numeric(ref_alt_strs[1])
            alt_rd = as.numeric(ref_alt_strs[2])
            total_rd = as.numeric(gt_strs[[1]][3])
            ref_ab = ref_rd/total_rd
            alt_ab = alt_rd/total_rd

            if (!is.nan(alt_ab)) {
                if (alt_rd >=minReads || alt_ab >=abCutoff){
                    gt_list = c(gt_list,gt)
                    ref_ab_list = c(ref_ab_list,ref_ab)
                    alt_ab_list = c(alt_ab_list,alt_ab)
                    s_list = c(s_list,s_name)
                }
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

sink(paste(outDir,".log_file.txt",sep="/"))
geneListFile = '/nfs/research/dunham/samples/ddd/data/gene_list/id_mr_list_chrMap_prune.txt'

varData = fread(raw_var,sep="\t",header=T,stringsAsFactors=F,quote="")
cat("--Data Loaded",'\n')
geneChrom = fread(geneListFile,sep="\t",header=F,stringsAsFactors=F,quote="")

geneList = subset(geneChrom,geneChrom$V2==chrNum)
#nHetsList = c()
#nHomsList = c()
df.gene.path = data.frame()


for (i in 1:dim(geneList)[1]){
    gene_name = geneList[i,1]
    #gene_name = "LARGE"
    cat('--GeneName:',as.vector(gene_name$V1),'\t',i,'\n')
    gene_id.grep = gene_name #paste("",geneList[i,1],"$",sep="")
    gene_id.app  = paste("^",geneList[i,1],"$",sep="")
    #gene_id.app = "^ABCA4$"

    df.gene.grep = subset(varData,grepl(gene_id.grep,varData$ANNO_SYMBOL))
  
    if (dim(df.gene.grep)[1] !=0) {

        #Extract Gene-HGNC related variants
        df.gene = df.gene.grep[apply(data.frame(df.gene.grep$ANNO_SYMBOL),
                               1,function(x){
                                   y = any(grepl(gene_id.app,strsplit(x,"[|]")[[1]]))
                                   return(y)
                               }
                               ),]
        #df.gene.path = df.gene

        if (dim(df.gene)[1] !=0) {
            #Extract Pathogenic and likely_pathogenic variants

            df.gene.vep.path = df.gene[apply(data.frame(df.gene[,c("ANNO_SYMBOL","ANNO_CLIN_SIG")]),1,function(x){
                                                        gene_list = x[1]
                                                        clin_sig_list = x[2]
                                                        gene_index = grep(gene_id.app,strsplit(gene_list,"[|]")[[1]])
                                                        clin_sig_strs = strsplit(clin_sig_list,"[|]")[[1]][gene_index]
                                                        clin_sig_strs_and = unlist(apply(data.frame(clin_sig_strs),
                                                                                  1,function(x){
                                                                                    y = strsplit(x,"[&]")[[1]]
                                                                                    }
                                                                                ))
                                                        y = any(grepl('pathogenic$|likely_pathogenic$',clin_sig_strs_and))
                                                        return(y)
                                                        }
                                            ),]
            df.gene.clinvar.path.1 = df.gene[apply(data.frame(df.gene[,c("ANNO_SYMBOL","CLNSIG")]),1,function(x){
                                                        gene_list = x[1]
                                                        clin_sig_list = x[2]
                                                        gene_index = grep(gene_id.app,strsplit(gene_list,"[|]")[[1]])
                                                        clin_sig_strs = strsplit(clin_sig_list,"[|]")[[1]][gene_index]
                                                        clin_sig_strs_and = unlist(apply(data.frame(clin_sig_strs),
                                                                                  1,function(x){
                                                                                    y = strsplit(x,"[/]")[[1]]
                                                                                    }
                                                                                 ))
                                                        y = any(grepl('Pathogenic|Likely_pathogenic',clin_sig_strs_and))
                                                        return(y)
                                                        }
                                            ),]

            if (chrNum !="Y") {
                df.gene.clinvar.path.2 = df.gene[apply(data.frame(df.gene[,c("ANNO_SYMBOL","CLNSIGCONF")]),1,function(x){
                                                        gene_list = x[1]
                                                        clin_sig_list = x[2]
                                                        gene_index = grep(gene_id.app,strsplit(gene_list,"[|]")[[1]])
                                                        clin_sig_strs = strsplit(clin_sig_list,"[|]")[[1]][gene_index]
                                                        clin_sig_strs_and = unlist(apply(data.frame(clin_sig_strs),
                                                                                  1,function(x){
                                                                                    y = strsplit(x,"[,]")[[1]]
                                                                                    }
                                                                                 ))
                                                        y = any(grepl('Pathogenic|Likely_pathogenic',clin_sig_strs_and))
                                                        return(y)
                                                        }
                                            ),]
            }else{
                df.gene.clinvar.path.2 = data.frame()
            }


            df.gene.path = unique(rbind(df.gene.vep.path,df.gene.clinvar.path.1,df.gene.clinvar.path.2))
            #df.gene.path = df.gene

            #Extract CADD >=25.0 variants
            #if (dim(df.gene.path)[1] !=0) {
            #    df.gene.cadd = df.gene.path[apply(data.frame(df.gene.path[,c("ANNO_SYMBOL","ANNO_CADD_PHRED")]),
            #                                    1,function(x){
            #                                        gene_list = x[1]
            #                                        cadd_list = x[2]
            #                                        gene_index = grep(gene_id.app,strsplit(gene_list,"[|]")[[1]])
            #                                        cadd_gene_strs = strsplit(cadd_list,"[|]")[[1]][gene_index]
            #                                        y = any(as.numeric(cadd_gene_strs) >=25.0)
            #                                        #y = any(cadd_gene_strs >=25.0)
            #                                        return(y)
            #                                    }
            #                         )
            #                   ,]
            #    df.gene.path = df.gene.cadd
            #}

        }else{
            df.gene.path = data.frame()
            df.gene.cadd = data.frame()
        }

        nHetsList = c()
        nHomsRefList = c()
        nHomsAltList = c()
        sHet_names_raw = c()
        sHet_names_RDAB = c()
        het_gt_raw = c()
        het_gt_RDAB = c()

        sample_index = which(colnames(df.gene.path)=="FORMAT")+1
        df.cols = dim(df.gene.path)[2]
        df.gene.path = as.data.frame(df.gene.path)
        #df.gene.path = as.data.frame(df.gene.cadd)

        if (dim(df.gene.path)[1]!=0){
            for (j in 1:dim(df.gene.path)[1]){
                
                sampleGT = df.gene.path[j,c(sample_index:df.cols)]

                tmp_nHets = apply(data.frame(sampleGT),1,function(x) grepl('^0/1:|^1/0:|^1:',x,perl=TRUE))

                sHet_names_all = unique(paste0(noquote(colnames(sampleGT)[tmp_nHets]),collapse=";"))
                het_gt_all = unique(paste0(noquote(sampleGT[,tmp_nHets]),collapse=";"))

                hetSampleGT = as.data.frame(sampleGT[,tmp_nHets])
                if (dim(hetSampleGT)[2]==1){
                    colnames(hetSampleGT) = sHet_names_all
                }else {
                    hetSampleGT = hetSampleGT
                }

                filterResult= filterRD_AB(hetSampleGT,5,0.15)
                s_het_name = filterResult[1]
                het_gt = filterResult[2]
                ref_ab_list = filterResult[3]
                alt_ab_list = filterResult [4] #filterRD_AB(sampleGT[,tmp_nHets],minReads=5,abCutoff=0.15)

                #nHets = length(which(tmp_nHets==TRUE))
                nHets_str = strsplit(het_gt,";")[[1]]
                
                if (length(nHets_str)==1 && nHets_str=="NA"){
                    nHets = 0
                }else{
                    nHets = length(nHets_str)
                }
                #Homzygous Ref
                tmp_nHomRef = apply(data.frame(sampleGT),1,function(x) grepl('^0/0:|^0:',x,perl=TRUE))
                nHomRef = length(which(tmp_nHomRef==TRUE))
                
                tmp_nHomAlt = apply(data.frame(sampleGT),1,function(x) grepl('^1/1:',x,perl=TRUE))
                nHomAlt = length(which(tmp_nHomAlt==TRUE))
                
                nHetsList = c(nHetsList,nHets)
                nHomsRefList = c(nHomsRefList,nHomRef)
                nHomsAltList = c(nHomsAltList,nHomAlt)
                sHet_names_raw = c(sHet_names_raw,sHet_names_all)
                sHet_names_RDAB = c(sHet_names_RDAB,s_het_name)
                het_gt_raw = c(het_gt_raw,het_gt_all)
                het_gt_RDAB = c(het_gt_RDAB,het_gt)
            }
            cat("--After for loop",'\n')
            #df.gene.path["nHets"] = nHetsList
            
            #df.gene.path["nHoms"] = nHomsList
            
            ####### Removes GT of all samples #########
            #df.gene.path.2 = cbind(df.gene.path[,c(1:(sample_index-1))],sHet_names_raw,sHet_names_RDAB,
            #                       het_gt_RDAB,nHetsList,nHomsRefList,nHomsAltList
            #                      )
            ####### Print GT of All Samples ########
            df.gene.path.2 = cbind(df.gene.path,sHet_names_raw,sHet_names_RDAB,
                                   het_gt_RDAB,nHetsList,nHomsRefList,nHomsAltList
                                  )
             
            ############################################
            
            df.gene.count = subset(df.gene.path.2,df.gene.path.2$nHetsList>0)
            out_file = paste(outDir,'/',gene_name,'_variants.txt',sep="")
            var_het_out_file = paste(outDir,'/',gene_name,'_variants_het.txt',sep="")
            count_file = paste(outDir,'/',gene_name,'_variants_count.txt',sep="")

            if (dim(df.gene.count)[1]>0) {
                var_cast_count = sum(df.gene.count$nHetsList)
                var_het_count = dim(df.gene.count)[1]
                count_file = paste(outDir,'/',gene_name,'_variants_count.txt',sep="")
                cat(as.character(gene_name),'\t',var_het_count,'\t',var_cast_count,'\n',file=count_file)
                write.table(df.gene.count,sep="\t",file=var_het_out_file,col.names=T,row.names=F,quote=F)
            }
            write.table(df.gene.path.2,sep="\t",file=out_file,col.names=T,row.names=F,quote=F)
        }
        else{
            cat('--No Pathogenic Variants found for:',gene_name$V1,'\n')
        }
    }else{
            cat('--Gene Donot exist:',gene_name$V1,'\n')
    }
}

#sink()
