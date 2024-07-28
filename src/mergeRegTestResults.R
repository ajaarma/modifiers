require('data.table')

out_path = '/hps/nobackup/dunham/ai-uk-can/analysis_geneset_regression'
#out_path = '/hps/nobackup/dunham/ai-uk-can/analysis_gene_burden'
res_path = '/hps/nobackup/dunham/ai-uk-can/analysis_geneset_regression/results'
#res_path = '/hps/nobackup/dunham/ai-uk-can/analysis_gene_burden/results'
cgt_test = c('sds','tmtA','tmtB','fit','dst','pmt','tdi','rtt','qualf')
var_class = c('frame','miss','syn')
filter_class = c('impact','urv')

c_method = "bonferroni"
#c_method = "fdr"
count = 1
cgt_list = list()


for(cgt in cgt_test){
    
    for(var_type in var_class){

        for (filt_ele in filter_class){
            impact_file =  paste(out_path,'/',cgt,'/',var_type,'/',var_type,'_',filt_ele,'_lasso_sort_sig.txt',sep="")
            #urv_file =  paste(out_path,'/',cgt,'/',var_type,'/',var_type,'_urv_lasso_sort_sig.txt',sep="")
            impData = read.table(impact_file,sep="\t",header=T)
            #urvData = read.table(urv_file,sep='\t',header=T)
        
            imp_nrow = dim(impData)[1]
            #urv_nrow = dim(urvData)[1]
        
            filt_imp_class = rep(filt_ele,imp_nrow) 
            var_imp_class  = rep(var_type,imp_nrow)
            cgt_imp_class = rep(cgt,imp_nrow)
            tmpDF = cbind(cgt_imp_class,var_imp_class,filt_imp_class,impData)
            count = count+1

            #p_value_adj = p.adjust(tmpDF$PValue,method="bonferroni")
            p_value_adj = p.adjust(tmpDF$PValue,method=c_method)
            col_p_val = which(colnames(tmpDF)=="PValue")
            tmpDF_1 = data.frame(tmpDF)[,c(1:col_p_val)]
            tmpDF_2 = data.frame(tmpDF[,((col_p_val+1):dim(tmpDF)[2])])

            tmpDF.adj = cbind(tmpDF_1,p_value_adj,tmpDF_2)
            colnames(tmpDF.adj)[11] = "SampleCount"
            cgt_list[[count]] = tmpDF.adj

        }
    }
}

if(c_method == 'bonferroni'){
    c_method = "bf"
}else if (c_method =="fdr"){
    c_method = c_method
}
mergeDF = rbindlist(cgt_list)
outFile.all = paste(res_path,"/",c_method,'/out_vars_all.txt',sep="")
write.table(mergeDF,file=outFile.all,sep='\t',row.names=F,col.names=T,quote=F)

#mergeDF.ncov = subset(mergeDF,!Predictors %in% c("age_lt50","age_gt65","AlcoholIntake","Gender","(Intercept)","ncars","rent","emp","nhh"))
mergeDF.ncov.1 = subset(mergeDF,!Predictors %in% c("age","age_sq","sex_age","sex_age_sq","AlcoholIntake","Gender","(Intercept)","ncars","rent","emp","nhh"))
mergeDF.ncov = subset(mergeDF.ncov.1,!grepl('f.22009',Predictors))
outFile.ncov = paste(res_path,'/',c_method,'/out_vars_ncov.txt',sep="")
write.table(mergeDF.ncov,file=outFile.ncov,sep='\t',row.names=F,col.names=T,quote=F)

mergeDF.0_05 = subset(mergeDF.ncov,p_value_adj <= 0.05)
outFile.0_05 = paste(res_path,'/',c_method,'/out_vars_0_05.txt',sep="")
write.table(mergeDF.0_05,file=outFile.0_05,sep='\t',row.names=F,col.names=T,quote=F)

mergeDF.0_01 = subset(mergeDF.ncov,p_value_adj <=0.01)
outFile.0_01 = paste(res_path,'/',c_method,'/out_vars_0_01.txt',sep="")
write.table(mergeDF.0_01,file=outFile.0_01,sep='\t',row.names=F,col.names=T,quote=F)

mergeDF.0_001 = subset(mergeDF.ncov,p_value_adj <0.001)
outFile.0_001 = paste(res_path,'/',c_method,'/out_vars_0_001.txt',sep="")
write.table(mergeDF.0_001,file=outFile.0_01,sep='\t',row.names=F,col.names=T,quote=F)

#a1 = mapTDI(mergeDF.0_001,'tdi')
