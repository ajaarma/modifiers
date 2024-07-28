################################################################################
# Description: Script to process all dnx-hpc variants (hpc_dnx_prior_ukbb) and
# divided them into three Case scenarios:
#   (a) Case1: Rare variants; gnomAD_AF <=1% CADD >=25.0; IHC >=6
#   (b) Case2: Ultra-rare; 0 <= gnomAD_AC <= 50; CADD >=25.0; IHC <=50
#   (c) Case3: Ultra-rare; 0 <= gnomAD_AC <= 50; CADD <25.0; IHC !=0
#
#
#
################################################################################

require('data.table')
require('R.utils')
require('dplyr')
require('Matrix')

createVarsDF <- function(geneFile,gene_id){
    #Subroutine to process all variants and add in a list
    varsData = fread(geneFile,sep='\t',header=T,stringsAsFactors=F,quote="")

    varsData.count = subset(varsData,SumCount >0)
    varsData.prior = subset(varsData.count,transPriorList !="NA")
    #varsData.prior = varsData

    #comb_samples_list = apply(data.frame(varsData.prior[,c('UniqSamples','hpc_samples_list')]),1,function(x){
    #                            uniq_samples = strsplit(x[1],';')[[1]]
    #                            hpc_samples = strsplit(x[2],';')[[1]]
    #                            comb_tmp = c(uniq_samples,hpc_samples)
    #                            comb_list = comb_tmp[which(comb_tmp!="NA")]
    #                            comb_ret = paste(unique(comb_list),sep=";")
    #                            return(comb_ret)
    #                        })
    #if (dim(varsData.prior)[1] ==1){
    #    comb_samples_list = paste(as.vector(comb_samples_list),sep=";")
    #}
    #return(comb_samples_list)

    if (dim(varsData.prior)[1] >0){
        gene_var_pos_list = apply(data.frame(varsData.prior[,'VAR_KEY']),1,function(x){
                                        tmp = x[1]
                                        tmp.2 = paste(gene_id,tmp,sep='_')
                                        return(tmp.2)
                                })

        cadd_phred_score = apply(data.frame(varsData.prior[,'ANNO_CADD_PHRED']),1,function(x){
                                        tmp = unique(strsplit(x,'[|]')[[1]])
                                        cadd_tmp = sort(tmp,decreasing=T)
                                        cadd_score = as.numeric(cadd_tmp[1])
                                        return(cadd_score)
                                })
       
        #varsData.2 = cbind(varsData.prior,comb_samples_list,cadd_phred_score,gene_name_list)
        varsData.2 = cbind(varsData.prior,cadd_phred_score,gene_var_pos_list)
    
        #colAnno = c('VAR_KEY','GNOMADe4_AC','GNOMADg4_AC','AC','AN','cadd_phred_score','comb_samples_list','SumCount')
        colAnno = c('VAR_KEY','GNOMADe4_AC','GNOMADg4_AC','GNOMADe4_AF','GNOMADg4_AF','AC','AN','AF',
                    'cadd_phred_score','CLNSIG','CLNSIGCONF','amPriorList','clnPriorList','revelPriorList',
                    'pliPriorList','consqPriorList','UniqSamples','hpc_samples_list','SumCount','gene_var_pos_list')

        #if (gene_id != 'SRY'){
        #    varSub = varsData.2[,c('VAR_KEY','GNOMADe4_AC','GNOMADg4_AC','AC','AN',
        #                       'cadd_phred_score','CLNSIG','CLNSIGCONF','amPriorList','clnPriorList',
        #                       'UniqSamples','hpc_samples_list','SumCount','gene_var_pos_list')]
        #}else {
        #    varSub = varsData.2[,c('VAR_KEY','GNOMADe4_AC','AC','AN',
        #                       'cadd_phred_score','CLNSIG','CLNSIGCONF','amPriorList','clnPriorList',
        #                       'UniqSamples','hpc_samples_list','SumCount','gene_var_pos_list')]
        varsData.2.cols = colnames(varsData.2)

        #varSub = varsData.2[,varsData.2.cols[varsData.2.cols %in% colAnno]]
        varSub = select(varsData.2,varsData.2.cols[varsData.2.cols %in% colAnno])
        varSub$geneID = as.vector(sapply(varSub$gene_var_pos_list,function(x){tmp = strsplit(x,"_")[[1]][1];return(tmp)}))
        varSub.path = varSub


        ###### Commented out to filter against Pathogenecity ######
        #if ('CLNSIGCONF' %in% colnames(varSub)) {
        #    varSub.path = subset(varSub,grepl('Pathogenic|Likely_pathogenic',CLNSIG) |
        #                            grepl('Pathogenic|Likely_pathogenic',CLNSIGCONF) |
        #                            grepl('pathogenic|likely_pathogenic',amPriorList) |
        #                            grepl('pathogenic|likely_pathogenic',amPriorList)
        #                     )
        #}else{
        #    varSub.path = subset(varSub,grepl('Pathogenic|Likely_pathogenic',CLNSIG) |
        #                            grepl('pathogenic|likely_pathogenic',amPriorList) |
        #                           grepl('pathogenic|likely_pathogenic',amPriorList)
        #                        )
        #}
        
    }else{
        varSub.path = data.frame()
        cat('Zero rows data frame: ',gene_id,'\n')
    }
    return(varSub.path)
}

getVarList <- function(geneVarsFileList,hpc_dnx_dir){

    #geneVarsFileList = geneVarsFileList[335]
    varsList = list()
    #for (i in 1:length(geneVarsFileList)){
    for (i in 1:length(geneVarsFileList)){
        geneFile = paste(hpc_dnx_dir,geneVarsFileList[i],sep='/')
        strs = strsplit(geneVarsFileList[i],"_")[[1]]
        gene_id = strs[1]
        cat(i,'\t',gene_id,'\n')   
        #varsList[[i]] = list(createVarsDF(geneFile,gene_id),gene_id)
        outDF = createVarsDF(geneFile,gene_id)
        if(dim(outDF)[1] !=0) {
            varsList[[i]] = createVarsDF(geneFile,gene_id)
        }
    }
    varCombAll = rbindlist(varsList,fill=T)
    outList = list()
    outList[[1]] = varsList
    outList[[2]] = varCombAll

    return(outList)
}

splitVarsByCases <- function(varData,case_id){
    #Subroutine to extract variants by Cases (See description)

    if(case_id == 'case1'){
        varCase.1 = subset(varData,SumCount >50 & cadd_phred_score >=25)
        varCase.2 = subset(varCase.1,AF <=0.01)
        varCase = varCase.2

        #Filter by LoF and Revel score
        #varCase.lof = subset(varCase,revelPriorList >=0.70 | pliPriorList >=0.70)
        varCase.lof = subset(varCase, pliPriorList >=0.70)

    }else if(case_id == 'case2'){
        varCase.na = subset(varData,is.na(GNOMADe4_AC) | is.na(GNOMADg4_AC))
        varCase.1 = subset(varData, (0 <= GNOMADe4_AC & GNOMADe4_AC <= 50 ) & (0 <= GNOMADg4_AC & GNOMADg4_AC <= 50))
        varCase.1.na = unique(rbind(varCase.1,varCase.na))
        varCase.2 = subset(varCase.1.na, SumCount <=50 & cadd_phred_score >=25)
        #varCase = unique(rbind(varCase.na,varCase.2))
        varCase = varCase.2
        
        #Filter by pli score >=0.90
        #varCase.lof = subset(varCase,revelPriorList >=0.70 | pliPriorList >=0.70)
        varCase.lof = subset(varCase, pliPriorList >=0.90)

    }else if (case_id == 'case3'){
        varCase.na = subset(varData,is.na(GNOMADe4_AC) | is.na(GNOMADg4_AC))
        varCase.1 = subset(varData, (0 <= GNOMADe4_AC & GNOMADe4_AC <= 50 ) & (0 <= GNOMADg4_AC & GNOMADg4_AC <= 50))
        varCase.1.na = unique(rbind(varCase.1,varCase.na))
        varCase.2 = subset(varCase.1, SumCount <=50 & cadd_phred_score <25)
        #varCase = unique(rbind(varCase.3,varCase.na))
        varCase = varCase.2

        #Filter by LoF and Revel Score
        #varCase.lof = subset(varCase,revelPriorList >=0.70 | pliPriorList >=0.70)
        varCase.lof = subset(varCase, pliPriorList >=0.70)
    }
    caseList = list()
    caseList[[1]] = data.frame(varCase)
    caseList[[2]] = data.frame(varCase.lof)
    
    return(caseList)
}

splitVarsByImpact <- function(varsData,anal_type){
    #Subroutine to filter variants by Impact
    varsData = subset(varsData,AF <=1e-5 | is.na(AF)) 
    if(anal_type == 'frame'){
        cat('--Processing Frameshift variant: ',anal_type,'\n')
        var.1 = subset(varsData,consqPriorList == 'frameshift_variant')
        var.2 = subset(varsData,consqPriorList == 'frameshift_variant&splice_region_variant')
        var.3 = subset(varsData,consqPriorList == 'frameshift_variant&start_lost')
        var.4 = subset(varsData,consqPriorList == 'splice_acceptor_variant')
        var.5 = subset(varsData,consqPriorList == 'splice_donor_variant')
        var.6 = subset(varsData,consqPriorList == 'stop_gained')
        var.7 = subset(varsData,consqPriorList == 'stop_gained&frameshift_variant')
        var.8 = subset(varsData,consqPriorList == 'stop_gained&frameshift_variant&splice_region_variant')
        var.9 = subset(varsData,consqPriorList == 'start_lost')

        var.frame = unique(rbind(var.1,var.2,var.3,var.4,var.5,var.6,var.7,var.8,var.9))
        var.frame.lof = subset(var.frame,pliPriorList >=0.90)
        return(var.frame.lof)
    
    }else if(anal_type == 'miss'){
        cat('--Processing Missense variant: ',anal_type,'\n')
        var.1 = subset(varsData,consqPriorList == 'missense_variant')
        var.2 = subset(varsData,consqPriorList == 'missense_variant&splice_region_variant')

        var.miss = unique(rbind(var.1,var.2))
        var.miss.cadd.lof = subset(var.miss,(cadd_phred_score >=25 | revelPriorList >=0.5 )& pliPriorList >=0.90)
        return(var.miss.cadd.lof)

    }else if(anal_type == 'syn'){
        cat('--Processing Synonymous variant: ',anal_type,'\n')
        var.1 = subset(varsData,consqPriorList == 'synonymous_variant')
        var.2 = subset(varsData,consqPriorList == 'splice_region_variant&synonymous_variant')
        var.syn = unique(rbind(var.1,var.2))
        var.syn.lof = subset(var.syn,pliPriorList >=0.90)

        return(var.syn.lof)
    }
}

splitVarsByURV <- function(varsData,anal_type){
    #Subroutine to filter variants by Impact
    varsData = subset(varsData,AF <=1e-5 | is.na(AF)) 
    
    if(anal_type == 'frame'){
        cat('--Processing Frameshift variant: ',anal_type,'\n')
        
        var.1 = subset(varsData,consqPriorList == 'frameshift_variant')
        var.2 = subset(varsData,consqPriorList == 'frameshift_variant&splice_region_variant')
        var.3 = subset(varsData,consqPriorList == 'frameshift_variant&start_lost')
        var.4 = subset(varsData,consqPriorList == 'splice_acceptor_variant')
        var.5 = subset(varsData,consqPriorList == 'splice_donor_variant')
        var.6 = subset(varsData,consqPriorList == 'stop_gained')
        var.7 = subset(varsData,consqPriorList == 'stop_gained&frameshift_variant')
        var.8 = subset(varsData,consqPriorList == 'stop_gained&frameshift_variant&splice_region_variant')
        var.9 = subset(varsData,consqPriorList == 'start_lost')

        var.frame = unique(rbind(var.1,var.2,var.3,var.4,var.5,var.6,var.7,var.8,var.9))
        #var.frame.lof = subset(var.frame,pliPriorList < 0.90 & cadd_phred_score >=25 )
        var.frame.lof = subset(var.frame,pliPriorList < 0.90 )
        return(var.frame.lof)
    
    }else if(anal_type == 'miss'){
        cat('--Processing missense variant: ',anal_type,'\n')
        
        var.1 = subset(varsData,consqPriorList == 'missense_variant')
        var.2 = subset(varsData,consqPriorList == 'missense_variant&splice_region_variant')

        var.miss = unique(rbind(var.1,var.2))
        #var.miss.cadd.lof = subset(var.miss,cadd_phred_score >=25 & pliPriorList <0.90 & SumCount <=5)
        var.miss.cadd.lof = subset(var.miss,(cadd_phred_score >=25 | revelPriorList >=0.50) & pliPriorList <0.90)
        return(var.miss.cadd.lof)
    }else if(anal_type == 'syn'){
        
        cat('--Processing Synonymous variant: ',anal_type,'\n')
        var.1 = subset(varsData,consqPriorList == 'synonymous_variant')
        var.2 = subset(varsData,consqPriorList == 'splice_region_variant&synonymous_variant')
        var.syn = unique(rbind(var.1,var.2))
        var.syn.lof = subset(var.syn,pliPriorList <0.90 )
        return(var.syn.lof)
    }
}


saveVarsByCases <- function(anal_type){
    #Subroutine to save variants by Case scenarios

    hpc_dnx_dir = '/hps/nobackup/dunham/ai-uk-can/analysis/hpc_dnx_prior_ukbb'
    geneVarsFileList = list.files(hpc_dnx_dir,pattern='.txt')
    geneVarFileDF = data.frame(geneVarsFileList,as.vector(sapply(geneVarsFileList,function(x){tmp = strsplit(x,"_")[[1]][1];return(tmp)})))
    colnames(geneVarFileDF) = c("geneFile","GeneID")

    id_gene_mono = '/nfs/research/dunham/samples/ddd/data/id_genes/out/dd_pa_mono_uniq.txt'
    idMonoData = read.table(id_gene_mono,sep='\t')
    
    monoIDgeneFileDF = geneVarFileDF[geneVarFileDF$GeneID %in% idMonoData$V1,]
    monoIDgeneList = monoIDgeneFileDF$geneFile

    #load('/hps/nobackup/dunham/ai-uk-can/analysis/cases_ukbb/mat/varCombAll.df')

    #outList = getVarList(geneVarsFileList,hpc_dnx_dir)
    outList = getVarList(monoIDgeneList,hpc_dnx_dir)
    varCombAll = outList[[2]]
    if (anal_type == 'cases'){
        save(outList,file='/hps/nobackup/dunham/ai-uk-can/analysis/cases_ukbb/mat/outList.list')
        save(varCombAll,file='/hps/nobackup/dunham/ai-uk-can/analysis/cases_ukbb/mat/varCombAll.df')
    }else if (anal_type %in% v('impact','urv')){

        save(outList,file='/hps/nobackup/dunham/ai-uk-can/analysis/matrices/outList.list')
        save(varCombAll,file='/hps/nobackup/dunham/ai-uk-can/analysis/matrices/varCombAll.df')
    }

    ############ Case DF ############

    if(anal_type == 'cases'){
        case1_list = splitVarsByCases(varCombAll,'case1')
        case2_list = splitVarsByCases(varCombAll,'case2')
        case3_list = splitVarsByCases(varCombAll,'case3')
        case1_list.gen.df = case1_list[[1]]
        case1_list.lof.df = case1_list[[2]]
        save(case1_list.gen.df,file='/hps/nobackup/dunham/ai-uk-can/analysis/cases_ukbb/mat/case1_list.gen.df')
        save(case1_list.lof.df,file='/hps/nobackup/dunham/ai-uk-can/analysis/cases_ukbb/mat/case1_list.lof.df')

        case2_list.gen.df = case2_list[[1]]
        case2_list.lof.df = case2_list[[2]]
        save(case2_list.gen.df,file='/hps/nobackup/dunham/ai-uk-can/analysis/cases_ukbb/mat/case2_list.gen.df')
        save(case2_list.lof.df,file='/hps/nobackup/dunham/ai-uk-can/analysis/cases_ukbb/mat/case2_list.lof.df')

        case3_list.gen.df = case3_list[[1]]
        case3_list.lof.df = case3_list[[2]]
        save(case3_list.gen.df,file='/hps/nobackup/dunham/ai-uk-can/analysis/cases_ukbb/mat/case3_list.gen.df')
        save(case3_list.lof.df,file='/hps/nobackup/dunham/ai-uk-can/analysis/cases_ukbb/mat/case3_list.lof.df')

    }else if (anal_type == "impact") {
        frame.df = splitVarsByImpact(varCombAll,'frame')
        miss.df = splitVarsByImpact(varCombAll,'miss')
        syn.df = splitVarsByImpact(varCombAll,'syn')
    
        save(frame.df,file='/hps/nobackup/dunham/ai-uk-can/analysis/impact_ukbb/mat/frame.df')
        save(miss.df,file='/hps/nobackup/dunham/ai-uk-can/analysis/impact_ukbb/mat/miss.df')
        save(syn.df,file='/hps/nobackup/dunham/ai-uk-can/analysis/impact_ukbb/mat/syn.df')
    }else if (anal_type == "urv"){
        frame.df = splitVarsByURV(varCombAll,'frame')
        miss.df = splitVarsByURV(varCombAll,'miss')
        syn.df = splitVarsByURV(varCombAll,'syn')

        save(frame.df,file='/hps/nobackup/dunham/ai-uk-can/analysis/urv_ukbb/mat/frame.df')
        save(miss.df,file='/hps/nobackup/dunham/ai-uk-can/analysis/urv_ukbb/mat/miss.df')
        save(syn.df,file='/hps/nobackup/dunham/ai-uk-can/analysis/impact_ukbb/mat/syn.df')
    }

}

createVarMat <- function(case.df,case_id,case_str){

    #Subroutine to create variant-by-sample matrices
    # Example: outList = createVarMat(frame.df,'frame','impact')
    #logFile = paste('/hps/nobackup/dunham/ai-uk-can/analysis/cases_ukbb',case_id,

    #print(paste('/hps/nobackup/dunham/ai-uk-can/analysis/cases_ukbb/mat/',case_id,'_list.',case_str,'.df',sep=""))
    #load(paste('/hps/nobackup/dunham/ai-uk-can/analysis/cases_ukbb/mat/',case_id,'_list.',case_str,'.df',sep=""))
   
    if (grepl(case_str,"case")) {
        outPath = paste('/hps/nobackup/dunham/ai-uk-can/analysis/cases_ukbb',case_id,sep='/')
    }else if (grepl(case_str,"impact")){
        outPath = paste('/hps/nobackup/dunham/ai-uk-can/analysis/impact_ukbb',case_id,sep='/')
    }else if (grepl(case_str,"urv")){
        outPath = paste('/hps/nobackup/dunham/ai-uk-can/analysis/urv_ukbb',case_id,sep='/')
    }
   
    print(outPath)
    logFile = paste(outPath,'/.logFile_',case_str,'.log',sep="")

    outFile = paste(outPath,'/',case_id,'_',case_str,'_list',sep="")
    
    sink(logFile)

    nonIDsamples = '/nfs/research/dunham/samples/ukbb/data/450k/nonIDSamples.txt'
    nidDF = data.frame(fread(nonIDsamples,sep='\t',stringsAsFactors=F,quote=""))

    idSamples = '/nfs/research/dunham/samples/ukbb/data/450k/onlyIDSamples.txt'
    idDF = data.frame(fread(idSamples,sep='\t',stringsAsFactors=F,quote=""))

    case.df = data.frame(case.df)
    caseSampleList = c()
    for (i in 1:dim(case.df)[1]){
        dnx_s = strsplit(case.df[i,'UniqSamples'],';')[[1]]
        hpc_s = strsplit(case.df[i,'hpc_samples_list'],';')[[1]]
        dnx_hpc = unique(c(dnx_s,hpc_s))
        dnx_hpc.nna = dnx_hpc[dnx_hpc !="NA" & !is.na(dnx_hpc)]
        caseSampleList = c(caseSampleList,dnx_hpc.nna)
    
    }

    case.s.nnid = setdiff(caseSampleList,nidDF$V1)
    case.s.nid = setdiff(case.s.nnid,idDF$V1)

    case.s.mat.row = unique(c(nidDF$V1,case.s.nid,case.s.nnid))
    
    #mat = Matrix(nrow=dim(nidDF)[1],ncol=dim(case.df)[1],data=0,sparse=T)
    mat = Matrix(nrow=length(case.s.mat.row),ncol=dim(case.df)[1],data=0,sparse=T)

    #rownames(mat) = nidDF$V1
    rownames(mat) = case.s.mat.row
    #colnames(mat) = case.df$gene_var_pos_list
    colnames(mat) = case.df$geneID

    dnxNoMatchList = c()

    print(dim(case.df))
    
    for (i in 1:dim(case.df)[1]){
    #for (i in 1:100){
        message(i)
        dnx_s  = strsplit(case.df[i,'UniqSamples'],';')[[1]]
        hpc_s = strsplit(case.df[i,'hpc_samples_list'],';')[[1]]
        dnx_hpc = unique(c(dnx_s,hpc_s))
        dnx_hpc.na = dnx_hpc[dnx_hpc !="NA" & !is.na(dnx_hpc)]
        
        ind = which(rownames(mat) %in% dnx_hpc.na ==TRUE)
        ind.ni = which(dnx_hpc.na %in% rownames(mat)  ==FALSE)
        dnx_hpc_na.ni = dnx_hpc.na[ind.ni]
        
        dnxNoMatchList = c(dnxNoMatchList,dnx_hpc_na.ni)

        if (length(ind) !=0) {
            mat[rownames(mat) %in% dnx_hpc.na,i] = 1
        }else {
            cat('Matching 0 samples in DNX for variant: ',colnames(mat)[i],'\n')
        }
        if(length(dnx_hpc_na.ni) >0){
            cat('No match samples in DNX for variants: ',colnames(mat)[i],'\t',dnx_hpc_na.ni,'\n')
        }
    }
    dnxNoMatchList = unique(dnxNoMatchList)
    outList = list()
    outList[[1]] = mat
    outList[[2]] = dnxNoMatchList
    
    sink()
    save(outList,file=outFile)
    return(outList)
}

createGeneMat <- function(case.df,case_id,case_str){

    #Subroutine to create variant-by-sample matrices
    # Example: outList = createGeneMat(frame.df,'frame','impact')
    #logFile = paste('/hps/nobackup/dunham/ai-uk-can/analysis/cases_ukbb',case_id,

    #print(paste('/hps/nobackup/dunham/ai-uk-can/analysis/cases_ukbb/mat/',case_id,'_list.',case_str,'.df',sep=""))
    #load(paste('/hps/nobackup/dunham/ai-uk-can/analysis/cases_ukbb/mat/',case_id,'_list.',case_str,'.df',sep=""))
   
    if (grepl(case_str,"case")) {
        outPath = paste('/hps/nobackup/dunham/ai-uk-can/analysis/cases_ukbb',case_id,sep='/')
    }else if (grepl(case_str,"impact")){
        outPath = paste('/hps/nobackup/dunham/ai-uk-can/analysis/impact_ukbb',case_id,sep='/')
    }else if (grepl(case_str,"urv")){
        outPath = paste('/hps/nobackup/dunham/ai-uk-can/analysis/urv_ukbb',case_id,sep='/')
    }
   
    print(outPath)
    logFile = paste(outPath,'/.logFile_',case_str,'.txt',sep="")

    outFile = paste(outPath,'/',case_id,'_',case_str,'_list',sep="")
    
    sink(logFile)

    nonIDsamples = '/nfs/research/dunham/samples/ukbb/data/450k/nonIDSamples.txt'
    nidDF = data.frame(fread(nonIDsamples,sep='\t',stringsAsFactors=F,quote=""))

    idSamples = '/nfs/research/dunham/samples/ukbb/data/450k/onlyIDSamples.txt'
    idDF = data.frame(fread(idSamples,sep='\t',stringsAsFactors=F,quote=""))

    case.df = data.frame(case.df)
    caseSampleList = c()
    for (i in 1:dim(case.df)[1]){
        dnx_s = strsplit(case.df[i,'UniqSamples'],';')[[1]]
        hpc_s = strsplit(case.df[i,'hpc_samples_list'],';')[[1]]
        dnx_hpc = unique(c(dnx_s,hpc_s))
        dnx_hpc.nna = dnx_hpc[dnx_hpc !="NA" & !is.na(dnx_hpc)]
        caseSampleList = c(caseSampleList,dnx_hpc.nna)
    
    }

    case.s.nnid = setdiff(caseSampleList,nidDF$V1)
    case.s.nid = setdiff(case.s.nnid,idDF$V1)

    case.s.mat.row = unique(c(nidDF$V1,case.s.nid,case.s.nnid))
    
    mat = Matrix(nrow=length(case.s.mat.row),ncol=length(unique(case.df$geneID)),data=0,sparse=T)

    rownames(mat) = case.s.mat.row
    colnames(mat) = unique(case.df$geneID)

    dnxNoMatchList = c()

    print(dim(case.df))
    
    #for (i in 1:dim(case.df)[1]){
    for (gene_id in unique(case.df$geneID)) {
        var.df = subset(case.df, geneID==gene_id)
        message(gene_id)
        dnx_hpc_list = c()
        for (i in 1:dim(var.df)[1]){

            dnx_s  = strsplit(var.df[i,'UniqSamples'],';')[[1]]
            hpc_s = strsplit(var.df[i,'hpc_samples_list'],';')[[1]]
            dnx_hpc_list = c(dnx_hpc_list,unique(c(dnx_s,hpc_s)))
        }
        dnx_hpc_list = unique(dnx_hpc_list)
        dnx_hpc.na = dnx_hpc_list[dnx_hpc_list !="NA" & !is.na(dnx_hpc_list)]
        
        ind = which(rownames(mat) %in% dnx_hpc.na ==TRUE)
        ind.ni = which(dnx_hpc.na %in% rownames(mat)  ==FALSE)
        dnx_hpc_na.ni = dnx_hpc.na[ind.ni]
        
        dnxNoMatchList = c(dnxNoMatchList,dnx_hpc_na.ni)
        
        j = which(colnames(mat)==gene_id)
        if (length(ind) !=0) {
            mat[rownames(mat) %in% dnx_hpc.na,j] = 1
        }else {
            cat('Matching 0 samples in DNX for variant: ',colnames(mat)[j],'\n')
        }
        if(length(dnx_hpc_na.ni) >0){
            cat('No match samples in DNX for variants: ',colnames(mat)[j],'\t',dnx_hpc_na.ni,'\n')
        }
    }
    dnxNoMatchList = unique(dnxNoMatchList)
    outList = list()
    outList[[1]] = mat
    outList[[2]] = dnxNoMatchList
    
    sink()
    save(outList,file=outFile)
    return(outList)
}

getSampleList <- function(case.df){

    sampleList = c()

    for (i in 1:dim(case.df)[1]){
        uniq_s = strsplit(case.df[i,'UniqSamples'],";")[[1]]
        hpc_s = strsplit(case.df[i,'hpc_samples_list'],";")[[1]]
        dnx_hpc = unique(c(uniq_s,hpc_s))
        sampleList = c(sampleList,dnx_hpc)
    }
    sampleList = unique(sampleList)

    return(sampleList)
}

getUkbMat <- function() {


    #Subroutine to create the UKB Matrix and reformat it's columns
    #ukbTab = '/nfs/research/dunham/samples/ukbb/data/cg_test/samples_cg_test_672066.tab'
    ukbTab = '/nfs/research/dunham/samples/ukbb/data/cg_test/samples_cg_test_672066_v2.tab'
    colsDesc = '/nfs/research/dunham/samples/ukbb/data/cg_test/col_ind_290523_desc.txt'

    ukbDF = fread(ukbTab,sep='\t',header=T,stringsAsFactors=F,quote="")
    colsDF = read.table(colsDesc,sep='\t',header=F)

    colnames(ukbDF) = colsDF$V1
    save(ukbDF,file='/hps/nobackup/dunham/ai-uk-can/analysis/matrices/ukbTab.df')

    gmf = '/nfs/research/dunham/samples/ukbb/data/cg_test/samples_cg_test_v3.txt'
    gmfData.1 = read.table(gmf,sep='\t',header=T)
    gmfData = gmfData.1[,c(1,3)]
    colnames(gmfData) = c('f.eid_samples','f.189.0.0_tdi')
    ukbDF.2 = dplyr::left_join(ukbDF,gmfData,by = join_by(f.eid_samples == f.eid_samples))
    ukbDF = ukbDF.2
    save(ukbDF,file='/hps/nobackup/dunham/ai-uk-can/analysis/matrices/ukb.cgt.gmf')


    return(ukbDF)
    
}

execURVstep <- function(anal_type){

    #load('/hps/nobackup/dunham/ai-uk-can/analysis/impact_ukbb/mat/varCombAll.df')
    load('/hps/nobackup/dunham/ai-uk-can/analysis/matrices/varCombAll.df')
    if(anal_type == "frame"){
        frame.df = splitVarsByURV(varCombAll,anal_type)
        save(frame.df,file="/hps/nobackup/dunham/ai-uk-can/analysis/urv_ukbb/mat/frame.df")
        outList = createGeneMat(frame.df,'frame','urv')
        return(frame.df)
    }else if(anal_type == 'miss'){
        miss.df = splitVarsByURV(varCombAll,anal_type)
        save(miss.df,file="/hps/nobackup/dunham/ai-uk-can/analysis/urv_ukbb/mat/miss.df")
        outList = createGeneMat(miss.df,'miss','urv')
        return(miss.df)
    }else if(anal_type == 'syn'){
        syn.df = splitVarsByURV(varCombAll,anal_type)
        save(syn.df,file="/hps/nobackup/dunham/ai-uk-can/analysis/urv_ukbb/mat/syn.df")
        outList = createGeneMat(syn.df,'syn','urv')
        return(syn.df)
    }
}

execImpactStep <- function(anal_type){

    #load('/hps/nobackup/dunham/ai-uk-can/analysis/impact_ukbb/mat/varCombAll.df')
    load('/hps/nobackup/dunham/ai-uk-can/analysis/matrices/varCombAll.df')
    if(anal_type == "frame"){
        frame.df = splitVarsByImpact(varCombAll,anal_type)
        save(frame.df,file="/hps/nobackup/dunham/ai-uk-can/analysis/impact_ukbb/mat/frame.df")
        outList = createGeneMat(frame.df,'frame','impact')
        return(frame.df)
    
    }else if(anal_type == 'miss'){
        miss.df = splitVarsByImpact(varCombAll,anal_type)
        save(miss.df,file="/hps/nobackup/dunham/ai-uk-can/analysis/impact_ukbb/mat/miss.df")
        outList = createGeneMat(miss.df,'miss','impact')
        return(miss.df)
    
    }else if(anal_type == 'syn'){
        syn.df = splitVarsByImpact(varCombAll,anal_type)
        save(syn.df,file="/hps/nobackup/dunham/ai-uk-can/analysis/impact_ukbb/mat/syn.df")
        outList = createGeneMat(syn.df,'syn','impact')
        return(syn.df)
    }
}
