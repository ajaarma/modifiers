#.libPaths("/hps/software/users/dunham/R_lib/4.0.3/")
require('data.table')
require('ggplot2')
require('scales')

cgtGmfHeatMap <- function(run_type) {

    #outDir = "/nfs/research/dunham/samples/ukbb/data/results"
    if (run_type == 'id') {
        outDir = "/hps/nobackup/dunham/ukbb/id_lasso_comb"
        consqDir = "/hps/nobackup/dunham/ukbb/filterSNV/id_genes_ukbb_merge_all/"
    }else if (run_type =='int') {
        outDir = "/hps/nobackup/dunham/ukbb/int_lasso_comb"
        consqDir = "/hps/nobackup/dunham/ukbb/filterSNV/int_genes_ukbb_merge_all/"
    }
    #outDir = "/Users/ajaykumar/data/AI-CAN-UK/ukbb/results/id_genes_results"
    #outDir = "/Users/ajaykumar/data/AI-CAN-UK/ukbb/results/int_genes_results"
    #test_str = c("cgt_pmt_res","cgt_fit_res","cgt_dst_res","cgt_sds_res","cgt_tmtA_res","cgt_tmtB_res")
    test_str = c("pmt","fit","dst","sds","tmtA","tmtB","qualf","occup",'hhi','tdi')

    case_list = c("case1","case2","case3")

    #b = c(-21,-10,-5,0,5,10,21)
    #b = c(-21,-10,-3,0,3,10,21)
    b = c(-20,-10,-3,-0.5,0.5,3,10,20)
    #cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    #col_list = c("red", "magenta3", "blue4")
    #col_list = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
    #col_list = c("red4","red3","red2","orange","blue1","blue2","blueviolet")
    #col_list = c("red4","red3","coral2","coral1","coral","skyblue","blue","blue1","blue2","blue4")
    col_list = c("red4","red3","red2","red1","red","skyblue","blue","blue1","blue2","blue4")
    #col_list = c("red4","red3","coral2","skyblue","blue1","blue2","blue4")

    file_str = c('_vt_rc25_hc6','_vt_urcgte25_hc6','_vt_urclte25_hc6')
    #varCG_001 = data.frame()
    #varCG_02 = data.frame() #setNames(data.frame(matrix(ncol=4,nrow=0)),c('Variants','Beta','Significance','Test'))
    #varCG_05 = data.frame() #setNames(data.frame(matrix(ncol=4,nrow=0)),c('Variants','Beta','Significance','Test'))

    varCGList = list()

    for (case_id in case_list) {
        varCG_001 = data.frame()
        varCG_02 = data.frame() #setNames(data.frame(matrix(ncol=4,nrow=0)),c('Variants','Beta','Significance','Test'))
        varCG_05 = data.frame() #setNames(data.frame(matrix(ncol=4,nrow=0)),c('Variants','Beta','Significance','Test'))
 
        for (test_id in test_str){
            #for (file_id in file_str){
            
            case.df = read.table(paste0(consqDir,case_id,'_map_variants_consq.txt',collapse=""),sep='\t',head=T)
            case.ren = apply(data.frame(case.df$Variants),1,function(y){
                                    strs = strsplit(y,'[_]')[[1]]
                                    vars = paste0(strs[3],'_',strs[1],'_',strs[2],collapse="")
                                    return(vars)
                                    }
                            )
            mapCase = cbind(case.ren,case.df$Consequences)
            colnames(mapCase) = c('Variants','Consequences')

            if (case_id =="case1") {
            #file_str = c("merge_all_vt_rc25_hc6_sort_sig.txt")
                file_id = '_vt_rc25_hc6'
                file_id_str = c("Case1_lasso_sort_sig.txt")
            }else if (case_id =="case2"){
                #file_str = c("merge_all_vt_urcgte25_hc6_sort_sig.txt")
                file_id = '_vt_urcgte25_hc6'
                file_id_str = c("Case2_lasso_sort_sig.txt")
            }else if (case_id =="case3") {
                #file_str = c("merge_all_vt_urclte25_hc6_sort_sig.txt")
                file_id = '_vt_urclte25_hc6'
                file_id_str = c("Case3_lasso_sort_sig.txt")
            }
         
            #Data per GMF category for given Case scenario
            if(test_id %in% c('qualf','occup','hhi')) {
                test_id_abs = toupper(paste0('GMF - ',test_id,collapse=""))
                file_name = paste(outDir,test_id,file_id,file_id_str,sep='/')
                varData = fread(file_name,sep='\t',header=T,stringsAsFactor=F,quote="")
                colnames(varData)[1] = c('Variables')
                colnames(varData)[6] = c('FDRsignificance')
                colnames(varData)[2] = c('LassoEst')
                colnames(varData)[3] = c('GLMEst')
                colnames(varData)[4] = c('GLM-StdError')
                colnames(varData)[5] = c('TValue')

                varGene = subset(varData, !grepl('Male|Female',varData$Variables))

                varGene$Variables = apply(data.frame(varGene$Variables),1,function(x) {
                                                          strs = strsplit(x,'_')[[1]]
                                                          return (paste(strs[1],strs[2],strs[3],sep="_"))
                                                       }
                                        )
                #Assign 0 significance to all the classes. Created this variable
                #for heatmap plotting
                varGene$FDRsignificance <- 0
                varGene$LassoEst = as.numeric(varGene$LassoEst)
                print(varGene)
                
            }else {
                if (test_id == 'tdi'){
                    test_id_abs = toupper(paste0('GMF - ',test_id,collapse=""))
                }else {
                    test_id_abs = toupper(paste0('CT - ',test_id,collapse=""))
                }
                file_name = paste(outDir,test_id,file_id,file_id_str,sep='/')
                varData = fread(file_name,sep='\t',header=T,stringsAsFactor=F,quote="")
                colnames(varData)[1] = c('Variables')
                colnames(varData)[6] = c('FDRsignificance')
                colnames(varData)[2] = c('LassoEst')
                colnames(varData)[3] = c('GLMEst')
                #varGene = subset(varData, !grepl('Intercept|Male|Female|Gender|Alcohol',varData$Variables))
                varGene = subset(varData, !grepl('Male|Female',varData$Variables))
                varGene$Variables = apply(data.frame(varGene$Variables),1,function(x) {
                                                          strs = strsplit(x,'_')[[1]]
                                                          return (paste(strs[1],strs[2],strs[3],sep="_"))
                                                       }
                                        )
            }
            print('outside')
            #test_id_abs = toupper(paste0(''test_id)
            
            varGenePrune001 = subset(varGene,varGene$FDRsignificance <0.001)
            varGenePrune02 = subset(varGene,varGene$FDRsignificance <0.02)
            varGenePrune05 = subset(varGene,varGene$FDRsignificance <0.05)
            varGenePrune001$Test = rep(test_id_abs,dim(varGenePrune001)[1])
            varGenePrune02$Test = rep(test_id_abs,dim(varGenePrune02)[1])
            varGenePrune05$Test = rep(test_id_abs,dim(varGenePrune05)[1])
      
            varCG_001 = rbind(varCG_001,varGenePrune001)
            varCG_02 = rbind(varCG_02,varGenePrune02)
            varCG_05 = rbind(varCG_05,varGenePrune05)
        }


        #print(varCG_02[1:5,])
        
        #Case 1
        if (case_id=="case1"){

            varCG_consq = getConsqMap(varCG_02,mapCase)
            varCG_consq_out_file = paste0(outDir,'/pathwayER/',case_id,'_var.txt',collapse="")
            write.table(varCG_consq,file=varCG_consq_out_file,sep='\t',col.names=T,row.names=F,quote=F)

            varCGList[[1]] = varCG_consq
            
            #(plot1 = ggplot(varCG_02,aes(y=Variables,x=Test,fill=Estimate))+
            (plot1 = ggplot(varCG_02,aes(y=Variables,x=Test,fill=LassoEst))+
            #(plot1 = ggplot(varCG_02,aes(y=Variables,x=Test,fill=GLMEst))+
                    geom_tile()+
                   labs(x="\n Cognitive Tests (CTs) and General Measures of Functioning (GMFs) ", y= "Predictors\n",
                   title = "Lasso estimates of predictors across CTs and GMFs (Case 1) \n",
                   fill = "Coefficient \n Estimate \n")+
                theme(plot.title = element_text(face="bold",hjust=0.5,size=33),
                    axis.title.x = element_text(face="bold",colour="black",size=26),
                    axis.title.y = element_text(face="bold",colour="black",size=26),
                    legend.title = element_text(face = "bold",colour="magenta",size=16),
                    legend.text = element_text(face = "bold",colour="magenta",size=12),
                    legend.key.height = unit(3,'cm'),
                    legend.key.width = unit(2,'cm'),
                    axis.text.y = element_text(face = "bold",colour="brown",size=16),
                    axis.text.x = element_text(face = "bold",colour="brown",size=16,angle=-45))+
                    geom_text(aes(label=round(LassoEst,4)),color="white",size=6)+
                    #geom_text(aes(label=round(GLMEst,4)),color="white",size=6)+
                    scale_fill_gradientn(limits=c(min(b),max(b)),
                                  colours=col_list,
                                  breaks=b,labels=format(b)
                                  )
                )
            out_file = paste(outDir,"Case_1_plot.png",sep="/")
            ggsave(plot1,file=out_file,width=22,height=22,dpi=300,device='png')
            #dev.off()
            
        }else if (case_id =="case2"){
            #Case 2
            varCG_consq = getConsqMap(varCG_001,mapCase)
            varCG_consq_out_file = paste0(outDir,'/pathwayER/',case_id,'_var.txt',collapse="")
            write.table(varCG_consq,file=varCG_consq_out_file,sep='\t',col.names=T,row.names=F,quote=F)

            varCGList[[2]] = varCG_consq
          #(plot1 = ggplot(varCG_001,aes(y=Variables,x=Test,fill=Estimate))+
            (plot1 = ggplot(varCG_001,aes(y=Variables,x=Test,fill=LassoEst))+
            #(plot1 = ggplot(varCG_02,aes(y=Variables,x=Test,fill=GLMEst))+
            #(plot1 = ggplot(varCG_001,aes(y=Variables,x=Test,fill=GLMEst))+
                    geom_tile()+
                    labs(x="\n Cognitive Tests (CTs) and General Measures of Functioning (GMFs) ", y= "Predictors\n",
                    title = "Predictors significant across CTs and GMFs (Case 2) \n",
                    fill = "Coefficient \n Estimate \n")+
                    theme(plot.title = element_text(face="bold",hjust=0.5,size=30),
                        axis.title.x = element_text(face="bold",colour="black",size=26),
                        axis.title.y = element_text(face="bold",colour="black",size=26),
                        legend.title = element_text(face = "bold",colour="magenta",size=16),
                        legend.text = element_text(face = "bold",colour="magenta",size=14),
                        legend.key.height = unit(2,'cm'),
                        legend.key.width = unit(2,'cm'),
                        axis.text.y = element_text(face = "bold",colour="brown",size=16),
                        axis.text.x = element_text(face = "bold",colour="brown",size=16,angle=-45))+
              
                        geom_text(aes(label=round(LassoEst,4)),color="white",size=6)+
                        #geom_text(aes(label=round(GLMEst,4)),color="white",size=4)+
                        scale_fill_gradientn(limits=c(min(b),max(b)),
                                     colours=col_list,
                                     breaks=b,labels=format(b)
                                    )
            )
            out_file = paste(outDir,"Case_2_plot.png",sep="/")
            ggsave(plot1,file=out_file,width=22,height=22,dpi=300,device='png')
            #dev.off()
          
        }else if (case_id=="case3"){
        #Case 3
            varCG_consq = getConsqMap(varCG_001,mapCase)
            varCG_consq_out_file = paste0(outDir,'/pathwayER/',case_id,'_var.txt',collapse="")
            write.table(varCG_consq,file=varCG_consq_out_file,sep='\t',col.names=T,row.names=F,quote=F)

            varCGList[[3]] = varCG_consq

            #(plot1 = ggplot(varCG_02,aes(y=Variables,x=Test,fill=LassoEst))+
            #(plot1 = ggplot(varCG_02,aes(y=Variables,x=Test,fill=GLMEst))+
            #(plot1 = ggplot(varCG_001,aes(y=Variables,x=Test,fill=GLMEst))+
            (plot1 = ggplot(varCG_001,aes(y=Variables,x=Test,fill=LassoEst))+
                geom_tile()+
                labs(x="\n Cognitive Tests (CTs) and General Measures of Functioning (GMFs) ", y= "Predictors\n",
                    title = "Predictors significant across CTs and GMFs (Case 3) \n",
                    fill = "Coefficient \n Estimate \n")+
                theme(plot.title = element_text(face="bold",hjust=0.5,size=30),
                    axis.title.x = element_text(face="bold",colour="black",size=26),
                    axis.title.y = element_text(face="bold",colour="black",size=26),
                    legend.key.height = unit(2,'cm'),
                    legend.key.width = unit(2,'cm'),
                    legend.title = element_text(face = "bold",colour="magenta",size=16),
                    legend.text = element_text(face = "bold",colour="magenta",size=14),
                    axis.text.y = element_text(face = "bold",colour="brown",size=16),
                    axis.text.x = element_text(face = "bold",colour="brown",size=16,angle=-45))+
                    geom_text(aes(label=round(LassoEst,4)),color="white",size=6)+
                    #geom_text(aes(label=round(GLMEst,4)),color="white",size=6)+
                    scale_fill_gradientn(limits=c(min(b),max(b)),
                                    colours=col_list,
                                    breaks=b,labels=format(b)
                                    )
            )
            out_file = paste(outDir,"Case_3_plot.png",sep="/")
            ggsave(plot1,file=out_file,width=22,height=22,dpi=300,device='png')
            #dev.off()
        }
    }
    return(varCGList)
}

getConsqMap <- function(varCG,mapCase){
    # Get filtered variants mapped to consequences (obtained from Supplementary
    # datas sheet with full list of variants before applying 
  
    mapCase = data.frame(mapCase)
    varCG.gene = subset(varCG,!grepl('Intercept|age|Alcohol|cp|nhh|hhi|ncars|Gender|rent|emp',Variables))
    print(head(mapCase))
    print(head(varCG.gene))
    out = mapCase[match(as.vector(varCG.gene$Variables),mapCase[,1]),]
    print(head(out))
    out.df = cbind(varCG.gene,out$Consequences)
    colnames(out.df)[8] = "Consequences"
    return(out.df)
}


corMatHeatMap <- function(run_type,anal_type,out_path="") {
    
   
    caseList = list()
    caseList[["_vt_rc25_hc6"]] = 'Case1'
    caseList[["_vt_urcgte25_hc6"]] = 'Case2'
    caseList[["_vt_urclte25_hc6"]] = 'Case3'

    #ggplot2 themes

    #My_theme = theme(
    #            axis.title.x = element_text(size = 16),
    #            axis.text.x = element_text(size = 10),
    #            axis.title.y = element_text(size = 10)
    #            )

    gList = list()
    if (anal_type=='id'){
        testPath = paste0('/hps/nobackup/dunham/ukbb/id_genes_cor/',run_type,collapse="")
    }else if(anal_type =='int'){
        testPath = paste0('/hps/nobackup/dunham/ukbb/int_genes_cor/',run_type,collapse="")
    }
    
    for(case_id in names(caseList)) {
        
        cat(run_type,'\t',case_id,'\n')
        objPath = paste0(testPath,'/',case_id,'/',caseList[[case_id]],'_mcor_pr_mat',collapse="")
        load(objPath)

        g = ggplot(data=mm.mcor,aes(x=Var1,y=Var2,fill=value))+geom_tile()+
            scale_fill_gradient('Cor-Coeff', limits=c(-1, 1), breaks = c(-1, 0.50, 0, 0.50, 1),  
                               low = "blue", high = "red")+
             ggtitle(paste0('Correlation b/w covariates in: ',run_type,'/',caseList[[case_id]],collapse=""))+
            theme(plot.title = element_text(hjust = 0.5,face="bold",size=18))+
            theme(axis.text.x=element_text(angle = -90, hjust = 0,size=8))+
            theme(axis.text.y=element_text(angle = 0,size=8))+
            labs(x="Covariates",y="Covariates")+
            theme(axis.title.x = element_text(face="bold",colour="black",size=14))+
            theme(axis.title.y = element_text(face="bold",colour="black",size=14))+
            theme(legend.title=element_text(face='bold',size=12))+
            theme(legend.text=element_text(face='bold',size=12))
        if(out_path!=''){
            out_plot = paste0(testPath,'/',case_id,'/',caseList[[case_id]],'_plot.png',collapse="")
            ggsave(out_plot,width=7,height=7,units="in")
        }
        gList[[caseList[[case_id]]]] = g
    }
    return(gList)
}

