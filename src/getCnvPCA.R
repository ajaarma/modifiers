#args=(commandArgs(TRUE))
require('irlba')
#-------------------------
# Description: Get CNV pca 
#
##########

set.seed(123)

#cmp = args[[1]]
cmp = 'dnx'

OUT_PATH = '/hps/nobackup/dunham/ai-uk-can/ukb200k/cnv_mat/clusters'
OUT_DIR = paste(OUT_PATH,cmp,sep='/')
system(paste0('mkdir -p ',OUT_DIR,sep=""))

matPath = '/hps/nobackup/dunham/ai-uk-can/ukb200k/cnv_mat'
server_type = c('hpc','dnx')

for (i in c(1:8)){
    matFile = paste(matPath,'/',cmp,'_cnv_mat_',i,sep="")
    load(matFile)
    message('-- cnv mat loaded iter: ',i)
    batchMat = batch.df[,-c(1:2)]
    batch.df= c()
    batchSeq = list()
    if (cmp == 'dnx') {
        message('-- Inside DNX seq loop')
        batchSeq[[1]] = c(1:5718)
        batchSeq[[2]] = c(5719:11437)
        batchSeq[[3]] = c(11438:17156)
        batchSeq[[4]] = c(17157:22875)
        batchSeq[[5]] = c(22875:dim(batchMat)[2])
    }else if(cmp == 'hpc'){
        batchSeq[[1]] = c(1:4600)
        batchSeq[[2]] = c(4601:9200)
        batchSeq[[3]] = c(9201:13800)
        batchSeq[[4]] = c(13801:18400)
        batchSeq[[5]] = c(18401:dim(batchMat)[2])
    }

    for(j in 1:5){
        col_ind = batchSeq[[j]]
        #mat = t(batchMat[,col_ind])
        st_time = Sys.time()
        mat.pca = prcomp_irlba(t(batchMat[,col_ind]),n=5,center=T,scale=F)
        end_time <- Sys.time()
        time_taken = end_time - st_time
        message('-- Time Spent: ',time_taken)
        prop.var.ind = which(summary(mat.pca)$importance[3,] <=0.85)
        mat = mat.pca$x[,prop.var.ind]
        k1 = kmeans(mat,centers=5)
        samples.df = data.frame(colnames(batchMat)[col_ind],k1$cluster)
        colnames(samples.df) = c('samples','cluster')
        
        outFile = paste(OUT_DIR,paste0('cnv_mat_',i,'_',j,'.txt',sep=""),sep='/')
        write.table(samples.df,sep='\t',file=outFile,row.names=F,col.names=F,quote=F)
    }
}
