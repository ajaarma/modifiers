# GenVisR package for lolliplot
require("trackViewer")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(rtracklayer)
require("data.table")

features <- GRanges("chr1", IRanges(c(1, 501, 1001), 
                                    width=c(120, 400, 405),
                                    names=paste0("block", 1:3)),
                    fill = c("#FF8833", "#51C6E6", "#DFA32D"),
                    height = c(0.02, 0.05, 0.08))

SNP <- c(10, 100, 105, 108, 400, 410, 420, 600, 700, 805, 840, 1400, 1402)
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, names=paste0("snp", SNP)),
                     color = sample.int(6, length(SNP), replace=TRUE),
                     score = sample.int(5, length(SNP), replace = TRUE))
lolliplot(sample.gr, features)


####### Plot Tracks based on gene name ############
extdata <- system.file("extdata", package="trackViewer", mustWork=TRUE)
filename = file.path(extdata, "fox2.bed")
optSty <- viewGene("HSPA8", filenames=filename, format="BED",
                   txdb=TxDb.Hsapiens.UCSC.hg19.knownGene,
                   org="org.Hs.eg.db")

##### Lolliplot plot for AGTR2 genes ############
varFile = "/Users/aak/Library/CloudStorage/GoogleDrive-aak@ebi.ac.uk/Other computers/My MacBook Pro/data/AI-CAN-UK/gdd_v3/lolplot/AGTR2_variants_het.txt"
annoFile = "/Users/aak/Library/CloudStorage/GoogleDrive-aak@ebi.ac.uk/Other computers/My MacBook Pro/data/AI-CAN-UK/gdd_v3/lolplot/Homo_sapiens.GRCh38.100.gtf.gz"

var.df = fread(varFile,sep="\t",header=T)
anno.df = fread(annoFile,skip=5)
colnames(anno.df) = c("chrom","model","type","start","end","uknown","strand","temp","rest")
chrNum = unique(var.df$CHROM)
anno.chr.exon.df = subset(anno.df,chrom ==strsplit(chrNum,'chr')[[1]][2] & type=='exon')

gene_name = paste0('"','AGTR2','"',collapse="")
anno.gene.chr.exon.df = anno.chr.exon.df[apply(as.data.frame(anno.chr.exon.df$rest),
                                                       1,function(x) {
                                                                  y = grepl(gene_name,x)
                                                                  return(y)
                                                      }
                                               ),]

dim(anno.gene.chr.exon.df)


