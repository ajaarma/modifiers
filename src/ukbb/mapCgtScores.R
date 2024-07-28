mapTDI <- function(cgt.df,cgt_type){
    #Subroutine to map f.eids to TDI scores
    
    load('/hps/nobackup/dunham/ai-uk-can/analysis/tdi/mat/cgt.tdi.mat')
    dnx_hpc_s_list = list()
    dnx_hpc_s_score = list()
    cgt.df = data.frame(cgt.df)

    for (i in 1:dim(cgt.df)[1]){
        print(i)
        dnx_s = strsplit(cgt.df[i,'UniqSamples'],';')[[1]]
        hpc_s = strsplit(cgt.df[i,'hpc_samples_list'],';')[[1]]
        dnx_hpc_s = unique(dnx_s,hpc_s)
        dnx_hpc_s_list[[i]] = paste(dnx_hpc_s,sep=";")
        score_list = c()
        
        for (ele in dnx_hpc_s){
            #print(ele)
            if(cgt_type == 'tdi') {
                #score = subset(cgt.df,f.eid_samples==ele)$tdi
                score = cgt.rm.df[cgt.rm.df$f.eid_samples==ele,'tdi']
            }else if(cgt_type=="sds"){
                score = cgt.rm.df[cgt.rm.df$f.eid_samples==ele,'sds']
                #score = subset(cgt.df,f.eid_samples==ele)$sds
            }else if(cgt_type=="fit"){
                score = cgt.rm.df[cgt.rm.df$f.eid_samples==ele,'fit']
                #score = subset(cgt.df,f.eid_samples==ele)$fit
            }else if(cgt_type == 'tmtA'){
                score = cgtrm..df[cgt.rm.df$f.eid_samples==ele,'tmtA']
                #score = subset(cgt.df,f.eid_samples==ele)$tmtA
            }else if(cgt_type == 'tmtB'){
                score = cgt.rm.df[cgt.rm.df$f.eid_samples==ele,'tmtB']
                #score = subset(cgt.df,f.eid_samples==ele)$tmtB
            }else if(cgt_type == 'dst'){
                score = cgt.rm.df[cgt.rm.df$f.eid_samples==ele,'dst']
                #score = subset(cgt.df,f.eid_samples==ele)$dst
            }else if(cgt_type == 'rtt'){
                score = cgt.rm.df[cgt.rm.df$f.eid_samples==ele,'rtt']
                #score = subset(cgt.df,f.eid_samples==ele)$rtt
            }else if(cgt_type == 'pmt'){
                score = cgt.rm.df[cgt.rm.df$f.eid_samples==ele,'pmt']
                #score = subset(cgt.df,f.eid_samples==ele)$pmt
            }else if(cgt_type == 'qualf'){
                score = cgt.rm.df[cgt.rm.df$f.eid_samples==ele,'qualf']
                #score = subset(cgt.df,f.eid_samples==ele)$qualf
            } 
            
            if (length(score)==0){
                score_list = c(score_list,'NA')
            }else{
                score_list = c(score_list,score)
            }
        }
        #score_list_vec = rbindlist(score_list)
        dnx_hpc_s_score[[i]] = paste(score_list,sep=';')
    }
    cgt.df$dnx_hpc_comb = rbindlist(dnx_hpc_s_list)
    cgt.df$score = dnx_hpc_s_score
    
    return(cgt.df)
    #return(dnx_hpc_s_score)

}
