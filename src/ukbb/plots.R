require('ggplot2')

plot_PRS <- function(run_id){
   
    
   #-------------------------------------------------
    # Subroutine to plot the effect of PRS on CGT score
    #--------------------------------------------------
   prs_res_file = '/Users/ajaykumar/data/AI-CAN-UK/ukbb/results/prs/edu_prs_rv_cgt.txt'
   glm.df = read.table(prs_res_file,sep='\t',header=T)
  
    #run_id = 'tdi'
    case.df = subset(glm.df,run_type==run_id)
    #case.sort = case.df[order(case.df$levels,decreasing=F),]
    case.df$xaxis = c(0.10,1.10,2.10,3.10,4.10,
                      0.35,1.35,2.35,3.35,4.35,
                      0.60,1.60,2.60,3.60,4.60,
                      0.85,1.85,2.85,3.85,4.85)
                      
    if (run_id == 'fit') {
        y_axis_title = c("Effect on FIT score (Beta)")
        legend_position = c(0.25,0.9)
    }else if (run_id == "nedu") {
        y_axis_title = c("Effect on NEDU years (Beta)")
        legend_position = c(0.25,0.9)     
    }else if (run_id == "tdi") {
      y_axis_title = c("Effect on TDI (Beta)")
      legend_position = c(0.3,0.9) 
    }
    p = ggplot(case.df,aes(x=xaxis,y=beta,group=case_id,color=case_id))+
        geom_point()+
        geom_errorbar(aes(ymin=beta-se,ymax=beta+se),width=0.2,position=position_dodge(0.05))+
        
        labs(title=element_blank(),x="PRS Quantiles",y=y_axis_title)+
        theme_classic()+
        scale_x_continuous(breaks=c(0.5,1.5,2.5,3.5,4.5), labels=c("0-20%","20-40%","40-60%","60-80%","80-100%"))+
        geom_vline(xintercept = c(1,2,3,4), linetype = "dashed", color="grey74")+
        scale_color_manual(values = rev(c('blue','violet','orange','red')),
                         labels = rev(c("PTV and missense non-carriers",
                                        "PTV non-carriers and carriers of missense variants",
                                        "PTV carriers and non-carriers of missense variants",
                                        "PTV carries and Missense carriers")))+
        theme(legend.position=legend_position)+
        theme(legend.title = element_blank())
    print(p)
    ggsave(p,file=paste0("/Users/ajaykumar/data/AI-CAN-UK/ukbb/results/prs/",run_id,"_prs.png"))
    #scale_y_continuous(limits = c(-0.5,0.8))
    
      
}

plot_burden_intolerance <- function() {
  
  #----------------------------------------
  # Subroutine to plot gene burden scores
  #----------------------------------------
  sum_burden_file = '/Users/ajaykumar/data/AI-CAN-UK/ukbb/results/sum_burden/out_vars_0_05.txt'
  gene_burden_file = '/Users/ajaykumar/data/AI-CAN-UK/ukbb/results/gene_burden/out_vars_0_05.txt'
  
  sum.df = read.table(sum_burden_file,sep='\t',header=T)
  gene.df = read.table(gene_burden_file,sep='\t',header=T)  
  
  for (cgt_type in factor(gene.df$cgt_imp_class)){
     
     if(cgt_type=="qualf"){
       gene.df$cgt_imp_class[gene.df$cgt_imp_class==cgt_type] = toupper("NEDU")
     }else if(cgt_type == "tmtA") {
       gene.df$cgt_imp_class[gene.df$cgt_imp_class==cgt_type] = toupper("TMT-A")
     }else if(cgt_type == "tmtB"){
       gene.df$cgt_imp_class[gene.df$cgt_imp_class==cgt_type] = toupper("TMT-B")
     }else{
       gene.df$cgt_imp_class[gene.df$cgt_imp_class==cgt_type] = toupper(cgt_type)
     }
  }

  cgt_label = c("NEDU","FIT","DST","SDS","PMT","RTT","TMT-A","TMT-B","TDI")
  
  data_summary <- function(x) {
    m <- median(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
  }
  
  #for (var_type in c("frame",""))
  gene.df.imp.frame = subset(gene.df,filt_imp_class=='impact'& var_imp_class=='frame')  
  gene.df.imp.miss = subset(gene.df,filt_imp_class=='impact'& var_imp_class=='miss')  
  gene.df.imp.syn = subset(gene.df,filt_imp_class=='impact'& var_imp_class=='syn')
  

  
  p = ggplot(subset(gene.df,filt_imp_class=='impact'& var_imp_class=='frame'),
             aes(Lasso.Est,y=factor(cgt_imp_class),fill="darkred"))+
    ggtitle("PTVs in Intolerant ID genes")+xlab("Beta Estimate (Lasso)")+ylab("Cognitive Tests")+
    #theme_classic()+
    geom_violin(scale="count",width=1.5)+geom_boxplot(width=0.15,color="black",fill="white")+
    scale_x_continuous(limits=c(-0.4,0.4))+ 
    scale_y_discrete(limits=cgt_label)+
    #stat_summary(fun.data=data_summary)+
    geom_vline(xintercept = c(0), linetype = "dashed", color="black")+
    #labs(x="LASSO/GLM (Beta) Estimate",y="Cognitive Tests",colour="black")+
    #scale_y_discrete(labels=cgt_label)+
    scale_fill_discrete(guide="none")+
    theme(axis.title.x = element_text(size = 14,vjust=-1.0,face="bold"), 
          axis.title.y = element_text(size = 14,hjust=0.5,face="bold"),
          text=element_text(size=12,face="bold"),
          plot.title=element_text(size=16,hjust=0.5,face="bold"))
    
  #print(p)
  ggsave(p,file="/Users/ajaykumar/data/AI-CAN-UK/ukbb/results/gene_burden/plots/PTV_intol_lasso.png")
    
  ####----- Missense Intolerant genes ######
  
  p = ggplot(subset(gene.df,filt_imp_class=='impact'& var_imp_class=='miss'),
             aes(Lasso.Est,y=factor(cgt_imp_class),fill="darkred"))+
    ggtitle("Missense in Intolerant ID genes")+xlab("Beta Estimate (LASSO)")+ylab("Cognitive Tests")+
    #theme_classic()+
    geom_violin(scale="count",width=1.5)+geom_boxplot(width=0.15,color="black",fill="white")+
    scale_x_continuous(limits=c(-0.4,0.4))+ 
    scale_y_discrete(limits=cgt_label)+
    #stat_summary(fun.data=data_summary)+
    geom_vline(xintercept = c(0), linetype = "dashed", color="black")+
    #labs(x="LASSO/GLM (Beta) Estimate",y="Cognitive Tests",colour="black")+
    #scale_y_discrete(labels=cgt_label)+
    scale_fill_discrete(guide="none")+
    theme(axis.title.x = element_text(size = 14,vjust=-1.0,face="bold"), 
          axis.title.y = element_text(size = 14,hjust=0.5,face="bold"),
          text=element_text(size=12,face="bold"),
          plot.title=element_text(size=16,hjust=0.5,face="bold"))
  #print(p)
  ggsave(p,file="/Users/ajaykumar/data/AI-CAN-UK/ukbb/results/gene_burden/plots/Miss_intol_lasso.png")
  ####----- Synonymous Intolerant genes ######
  
  p = ggplot(subset(gene.df,filt_imp_class=='impact'& var_imp_class=='syn'),
             aes(GLM.Est,y=factor(cgt_imp_class),fill="darkred"))+
    ggtitle("Synonymous variants in Intolerant ID genes")+xlab("Beta Estimate (LASSO)")+ylab("Cognitive Tests")+
    #theme_classic()+
    geom_violin(scale="count",width=1.5)+geom_boxplot(width=0.15,color="black",fill="white")+
    scale_x_continuous(limits=c(-0.4,0.4))+ 
    scale_y_discrete(limits=cgt_label)+
    #stat_summary(fun.data=data_summary)+
    geom_vline(xintercept = c(0), linetype = "dashed", color="black")+
    #labs(x="LASSO/GLM (Beta) Estimate",y="Cognitive Tests",colour="black")+
    #scale_y_discrete(labels=cgt_label)+
    scale_fill_discrete(guide="none")+
    theme(axis.title.x = element_text(size = 14,vjust=-1.0,face="bold"), 
          axis.title.y = element_text(size = 14,hjust=0.5,face="bold"),
          text=element_text(size=12,face="bold"),
          plot.title=element_text(size=16,hjust=0.5,face="bold"))
  print(p)
  ggsave(p,file="/Users/ajaykumar/data/AI-CAN-UK/ukbb/results/gene_burden/plots/Syn_intol_lasso.png")
}


plot_burden_tolerance <- function() {
  
  #----------------------------------------
  # Subroutine to plot gene burden scores
  #----------------------------------------
  sum_burden_file = '/Users/ajaykumar/data/AI-CAN-UK/ukbb/results/sum_burden/out_vars_0_05.txt'
  gene_burden_file = '/Users/ajaykumar/data/AI-CAN-UK/ukbb/results/gene_burden/out_vars_0_05.txt'
  
  sum.df = read.table(sum_burden_file,sep='\t',header=T)
  gene.df = read.table(gene_burden_file,sep='\t',header=T)  
  
  for (cgt_type in factor(gene.df$cgt_imp_class)){
    
    if(cgt_type=="qualf"){
      gene.df$cgt_imp_class[gene.df$cgt_imp_class==cgt_type] = toupper("NEDU")
    }else if(cgt_type == "tmtA") {
      gene.df$cgt_imp_class[gene.df$cgt_imp_class==cgt_type] = toupper("TMT-A")
    }else if(cgt_type == "tmtB"){
      gene.df$cgt_imp_class[gene.df$cgt_imp_class==cgt_type] = toupper("TMT-B")
    }else{
      gene.df$cgt_imp_class[gene.df$cgt_imp_class==cgt_type] = toupper(cgt_type)
    }
  }
  
  cgt_label = c("NEDU","FIT","DST","SDS","PMT","RTT","TMT-A","TMT-B","TDI")
  
  data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
  }

  p = ggplot(subset(gene.df,filt_imp_class=='urv'& var_imp_class=='frame'),
             aes(Lasso.Est,y=factor(cgt_imp_class),fill="darkred"))+
    ggtitle("PTVs in Tolerant ID genes")+xlab("Beta Estimate (LASSO)")+ylab("Cognitive Tests")+
    #theme_classic()+
    geom_violin(scale="count",width=1)+geom_boxplot(width=0.15,color="black",fill="white")+
    scale_x_continuous(limits=c(-0.4,0.4))+ 
    scale_y_discrete(limits=cgt_label)+
    #stat_summary(fun.data=data_summary)+
    geom_vline(xintercept = c(0), linetype = "dashed", color="black")+
    #labs(x="LASSO/GLM (Beta) Estimate",y="Cognitive Tests",colour="black")+
    #scale_y_discrete(labels=cgt_label)+
    scale_fill_discrete(guide="none")+
    theme(axis.title.x = element_text(size = 14,vjust=-1.0,face="bold"), 
          axis.title.y = element_text(size = 14,hjust=0.5,face="bold"),
          text=element_text(size=12,face="bold"),
          plot.title=element_text(size=16,hjust=0.5,face="bold"))
  
  print(p)
  ggsave(p,file="/Users/ajaykumar/data/AI-CAN-UK/ukbb/results/gene_burden/plots/2_PTV_tol_lasso.png")
  ####----- Missense Intolerant genes ######
  
  p = ggplot(subset(gene.df,filt_imp_class=='urv'& var_imp_class=='miss'),
             aes(Lasso.Est,y=factor(cgt_imp_class),fill="darkred"))+
    ggtitle("Missense in Tolerant ID genes")+xlab("Beta Estimate (LASSO)")+ylab("Cognitive Tests")+
    #theme_classic()+
    geom_violin(scale="count",width=0.5)+geom_boxplot(width=0.15,color="black",fill="white")+
    scale_x_continuous(limits=c(-0.4,0.4))+ 
    scale_y_discrete(limits=cgt_label)+
    #stat_summary(fun.data=data_summary)+
    geom_vline(xintercept = c(0), linetype = "dashed", color="black")+
    #labs(x="LASSO/GLM (Beta) Estimate",y="Cognitive Tests",colour="black")+
    #scale_y_discrete(labels=cgt_label)+
    scale_fill_discrete(guide="none")+
    theme(axis.title.x = element_text(size = 14,vjust=-1.0,face="bold"), 
          axis.title.y = element_text(size = 14,hjust=0.5,face="bold"),
          text=element_text(size=12,face="bold"),
          plot.title=element_text(size=16,hjust=0.5,face="bold"))
  
  print(p)
  ggsave(p,file="/Users/ajaykumar/data/AI-CAN-UK/ukbb/results/gene_burden/plots/2_Miss_tol_lasso.png")
  ####----- Synobymous Intolerant genes ######
  
  p = ggplot(subset(gene.df,filt_imp_class=='urv'& var_imp_class=='syn'),
             aes(Lasso.Est,y=factor(cgt_imp_class),fill="darkred"))+
    ggtitle("Synonymous variants in Tolerant ID genes")+xlab("Beta Estimate (LASSO)")+ylab("Cognitive Tests")+
    #theme_classic()+
    geom_violin(scale="count",width=1.5)+geom_boxplot(width=0.15,color="black",fill="white")+
    scale_x_continuous(limits=c(-0.4,0.4))+ 
    scale_y_discrete(limits=cgt_label)+
    #stat_summary(fun.data=data_summary)+
    geom_vline(xintercept = c(0), linetype = "dashed", color="black")+
    #labs(x="LASSO/GLM (Beta) Estimate",y="Cognitive Tests",colour="black")+
    #scale_y_discrete(labels=cgt_label)+
    scale_fill_discrete(guide="none")+
    theme(axis.title.x = element_text(size = 14,vjust=-1.0,face="bold"), 
          axis.title.y = element_text(size = 14,hjust=0.5,face="bold"),
          text=element_text(size=12,face="bold"),
          plot.title=element_text(size=16,hjust=0.5,face="bold"))
  print(p)
  ggsave(p,file="/Users/ajaykumar/data/AI-CAN-UK/ukbb/results/gene_burden/plots/2_Syn_tol_lasso.png")
}
