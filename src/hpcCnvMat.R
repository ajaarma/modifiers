args=(commandArgs(TRUE))
require('data.table')
require('dplyr')

list_num = args[[1]]
#list_num = 1

cnv_path = '/hps/nobackup/dunham/ai-uk-can/ukb200k/hpc_cnv'
test_sample = '/hps/nobackup/dunham/ai-uk-can/ukb200k/hpc_cnv/batch_01/SRC/1028673_500.tsv'
a1 = read.table(test_sample,header=T, comment.char="@")
a2 = data.frame(paste(a1$CONTIG,a1$START,a1$END,sep="_"),a1$COUNT)
colnames(a2) = c("CHR_POS","dummy")
batch.df = c()
batch.df = a2

batchList = list.dirs(cnv_path,recursive=F,full.names=F)
batchList_split = split(batchList,rep_len(1:8,length(batchList)))

count = 1
for (batch_id in batchList_split[[list_num]]){
    cat(batch_id,'\t',count,'\t',format(Sys.time(), "%a %b %d %X %Y"),'\n')
    sample_path = paste0(cnv_path,'/',batch_id,'/SRC/',sep="")
    #a1 = read_bulk(sample_path,extension=".tsv",fun=fread,sep="\t",header=T,stringsAsFactors=F)
    sampleList = list.files(sample_path)
    for (s_id in sampleList){
        s_name = strsplit(s_id,"[_]")[[1]][1]
        s_id_path = paste0(sample_path,s_id,sep="")
        a1 = read.table(s_id_path,header=T, comment.char="@")
        a2 = data.frame(paste(a1$CONTIG,a1$START,a1$END,sep="_"),a1$COUNT)
        colnames(a2) = c("CHR_POS",s_name)
        #batch.df = dplyr::left_join(batch.df,a2,by=join_by(CHR_POS==CHR_POS))
        batch.df = cbind(batch.df,a2[,2])
        colnames(batch.df)[dim(batch.df)[2]] = s_name
    }
    count = count+1
}

save(batch.df,file=paste("/hps/nobackup/dunham/ai-uk-can/ukb200k/cnv_mat/hpc_cnv_mat_",list_num,sep=""))
