args=(commandArgs(TRUE))
require('data.table')
require('dplyr')

list_num = args[[1]]
#list_num = 1

cnv_path = '/hps/nobackup/dunham/ai-uk-can/ukb200k/dx_cnv/280k_cnv'
#test_sample = '/hps/nobackup/dunham/ai-uk-can/ukb200k/hpc_cnv/batch_01/SRC/1028673_500.tsv'
test_sample = '/hps/nobackup/dunham/ai-uk-can/ukb200k/dx_cnv/280k_cnv/dnx_1/slice_00/5991126_23143_0_0.tsv.gz'
a1 = read.table(test_sample,header=T, comment.char="@")
a2 = data.frame(paste(a1$CONTIG,a1$START,a1$END,sep="_"),a1$COUNT)
colnames(a2) = c("CHR_POS","dummy")

batch.df = c()
batch.df = a2

#batchList = list.dirs(cnv_path,recursive=F,full.names=F)
#batchList_split = split(batchList,rep_len(1:8,length(batchList)))
load('/hps/nobackup/dunham/ai-uk-can/ukb200k/cnv_mat/dnx_batchList_split_list')

count = 1
for (s_id in batchList_split[[list_num]]){
    #cat(s_id,'\t',count,'\t',format(Sys.time(), "%a %b %d %X %Y"),'\n')
    cat(count,'\n',file=paste('/hps/nobackup/dunham/ai-uk-can/ukb200k/cnv_mat/.dnx_mat_log_',list_num,'.txt',sep=""))
    s_name = strsplit(s_id,"[/]|[_]")[[1]][5]
    s_id_path = paste(cnv_path,s_id,sep="/")
    a1 = read.table(s_id_path,header=T, comment.char="@")
    a2 = data.frame(paste(a1$CONTIG,a1$START,a1$END,sep="_"),a1$COUNT)
    colnames(a2) = c("CHR_POS",s_name)
    #batch.df = dplyr::left_join(batch.df,a2,by=join_by(CHR_POS==CHR_POS))
    batch.df = cbind(batch.df,a2[,2])
    colnames(batch.df)[dim(batch.df)[2]] = s_name
    count = count+1
}

save(batch.df,file=paste("/hps/nobackup/dunham/ai-uk-can/ukb200k/cnv_mat/dnx_cnv_mat_",list_num,sep=""))
