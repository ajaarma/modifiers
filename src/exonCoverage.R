loc4_0 = "/hps/software/users/dunham/R_lib/4.0.3"
require('data.table',lib.loc=loc4_0)
require("optparse",lib.loc=loc4_0)

exon_list = c()
mean_cov = c()
median_cov = c()

option_list = list(
                     make_option(c("-v", "--vars"), type="character", default=NULL,
                                 help="variants file name", metavar="character"),
                     make_option(c("-o", "--out"), type="character", default="out.filt.txt",
                                 help="output file name [default= %default]", metavar="character")
                     )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$vars) || is.null(opt$out)){
      print_help(opt_parser)
  stop("At least two arguments must be supplied (variants file and transcript lengths).n", call.=FALSE)
}

chr_dir = opt$vars
out_path = opt$out

zero_size_list = c()
#for (chr_dir in list.dirs(analysis_dir,recursive=F)){

getSummaryStat = function(chr_dir) {
    out_frame = data.frame()
    counter = 0
    #if (grepl("chrY",chr_dir)) {
        for(exon_iter in list.files(chr_dir)){
            exon_file = paste(chr_dir,exon_iter,sep='/')
            
            if(file.size(exon_file)==0){
                #cat(exon_iter,file.size(exon_file),"\n")
                zero_size_list = c(zero_size_list,exon_iter)
            }else {
                counter = counter+1
                a1 = c()
                cat(counter,exon_iter,file.size(exon_file),"\n")
                zero_size_list = c(zero_size_list,exon_iter)
                a1 = read.table(exon_file,sep='\t')
                
                a11.mean = apply(a1[,2:dim(a1)[2]],1, 
                                    function(x) {
                                                  x1 = as.numeric(x);
                                                  y = mean(x1[!is.na(x1)]);
                                                  return(y)
                                                }
                                )         
                a11.median = apply(a1[,2:dim(a1)[2]],1, 
                                   function(x) {
                                                x1 = as.numeric(x);
                                                 y = median(x1[!is.na(x1)]);
                                                return(y)
                                               }
                                  )
                a11.mean.2 = round(mean(a11.mean),4)
                a11.median.2 = median(a11.median)
                exon_list = c(exon_list,strsplit(exon_iter,".txt")[[1]][1])
                mean_cov = c(mean_cov,a11.mean.2)
                median_cov = c(median_cov,a11.median.2)

            }
        }
        out.frame = as.data.frame(cbind(exon_list,mean_cov,median_cov))
    #}
    return(out.frame)
}

mat.frame = getSummaryStat(chr_dir)
#mat.frame = as.data.frame(cbind(exon_list,mean_cov,median_cov))
write.table(mat.frame,file=paste(out_path,paste0(basename(chr_dir),'_mat_frame.txt'),sep='/'),
            sep='\t',row.names=F,col.names=T,quote=F)

