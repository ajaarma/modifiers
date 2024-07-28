require('moments')


getDist <- function(cgt.rm.df.raw,test_id){
    #Subroutine to get the distribution of CGT test scores

    if(test_id == 'sds'){
        plotPath = "/hps/nobackup/dunham/ai-uk-can/analysis/sds/plots"
        plotFile_1 = paste(plotPath,'/case1_gen_raw.png',sep='')
        plotFile_2 = paste(plotPath,'/case1_gen_mod.png',sep='')

        cgt.rm.df = subset(cgt.rm.df.raw, (sds >3 & sds <33))

        outFile_txt = paste(plotPath,'/sds_stats.txt',sep='')
        outCon <- file(outFile_txt)
        kurt_raw = kurtosis(cgt.rm.df.raw$sds,na.rm=T)
        skew_raw = skewness(cgt.rm.df.raw$sds,na.rm=T)

        kurt_mod = kurtosis(cgt.rm.df$sds,na.rm=T)
        skew_mod = skewness(cgt.rm.df$sds,na.rm=T)

        outFile_txt = paste(plotPath,'/sds_stats.txt',sep='')
        outCon <- file(outFile_txt)
        writeLines(c('Kurtosis of Raw data: ',kurt_raw,'\n'),outCon)
        writeLines(c('Skewness of Raw data: ',skew_raw,'\n'),outCon)
        writeLines(c('Kurtosis of mod data: ',kurt_mod,'\n'),outCon)
        writeLines(c('skewness of mod data: ',skew_mod,'\n'),outCon)
        close(outCon)

        png(plotFile_1)
        hist(cgt.rm.df.raw$sds,breaks=150,main="Histogram plot for SDS raw data",
             xlab="SDS score","Frequency",freq=T,col="blue")
        dev.off()
        png(plotFile_2)
        hist(cgt.rm.df$sds,breaks=150,main="Histogram plot for SDS modified data",
             xlab="SDS score",ylab="Frequency",freq=T,col="blue",xlim=range(0:50))
        dev.off()
    
    }else if(test_id == 'tmtA'){
        plotPath = "/hps/nobackup/dunham/ai-uk-can/analysis/tmta/plots"
        plotFile_1 = paste(plotPath,'/case1_gen_raw.png',sep='')
        plotFile_2 = paste(plotPath,'/case1_gen_mod.png',sep='')
        
        outFile_txt = paste(plotPath,'/sds_stats.txt',sep='')
        outCon <- file(outFile_txt)
        kurt_raw = kurtosis(cgt.rm.df.raw$sds,na.rm=T)
        skew_raw = skewness(cgt.rm.df.raw$sds,na.rm=T)

        kurt_mod = kurtosis(cgt.rm.df$sds,na.rm=T)
        skew_mod = skewness(cgt.rm.df$sds,na.rm=T)

        outFile_txt = paste(plotPath,'/sds_stats.txt',sep='')
        outCon <- file(outFile_txt)
        writeLines(c('Kurtosis of Raw data: ',kurt_raw,'\n'),outCon)
        writeLines(c('Skewness of Raw data: ',skew_raw,'\n'),outCon)
        writeLines(c('Kurtosis of mod data: ',kurt_mod,'\n'),outCon)
        writeLines(c('skewness of mod data: ',skew_mod,'\n'),outCon)
        close(outCon)

        png(plotFile_1)
        hist(cgt.rm.df.raw$tmtA,breaks=150,main="Histogram plot for TMT-A (log) data",
             xlab="TMT-A score","Frequency",freq=T,col="blue")
        dev.off()
         
    }else if (run_type =='tmtb'){

        Y.tmtB.log = cgtrm.df.raw$tmtB
        hist(Y.tmtB.log[!is.na(Y.tmtB.log)],breaks=150,main="Histogram plot for TMT-B (log) data",xlab="TMT-B score","Frequency",freq=T,col="blue")
    
    }else if (run_type == 'fit'){

        hist(cgt.rm.df.raw$fit,breaks=50,main="Histogram plot for FIT data",xlab="FIT score","Frequency",freq=T,col="blue")
    
    }else if (run_type == 'dst'){
        
        hist(cgt.rm.df.raw$dst,breaks=15,main="Histogram plot for DST data",xlab="DST score","Frequency",freq=T,col="blue")
    }else if(run_type == 'pmt'){
        par(mfrow=c(2,3))
        hist(as.numeric(Y.pmt.vec.all[,1]),breaks=100,main="PMT Score Round-1",xlab="Round - 1",ylab="Frequency",col="blue")
        hist(as.numeric(Y.pmt.vec.all[,1]),breaks=100,main="PMT Score Round-1",xlab="Round - 1",ylab="Frequency",col="blue")
        hist(as.numeric(Y.pmt.vec.all[,2]),breaks=100,main="PMT Score Round-2",xlab="Round-2",ylab="Frequency",col="blue")
        hist(as.numeric(Y.pmt.vec.all[,3]),breaks=100,main="PMT Score Round-3",xlab="Round-3",ylab="Frequency",col="blue")
        hist(as.numeric(Y.pmt.vec.all[,4]),breaks=100,main="PMT Score All-Round-Sum",xlab="All-Round-Sum",ylab="Frequency",col="blue")
        hist(as.numeric(Y.pmt.vec.all[,5]),breaks=100,main="PMT Score All-Round-Sum-LogT",xlab="All-Round-Sum-LogT",ylab="Frequency",col="blue")
        hist(as.numeric(pmt.cov.vars.nna.nz$pmt),breaks=100,main="PMT Score All-Round-Sum-LogT-nz",xlab="All-Round-Sum-LogT-nz",ylab="Frequency",col="blue")
    }



}
