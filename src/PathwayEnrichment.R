require("pathfindR")

bd_genes_ = read.table("/Users/aak/data/AI-CAN-UK/gdd_v3/ms_gene_var_count_sig.txt")

if(!"rWikiPathways" %in% installed.packages()){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("rWikiPathways", update = FALSE)
}
library(rWikiPathways)

load.libs <- c(
  "DOSE",
  "GO.db",
  "GSEABase",
  "org.Hs.eg.db",
  "clusterProfiler",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "RColorBrewer",
  "rWikiPathways",
  "RCy3")
options(install.packages.check.source = "no")
options(install.packages.compile.from.source = "never")
if (!require("pacman")) install.packages("pacman"); library(pacman)
p_load(load.libs, update = TRUE, character.only = TRUE)
status <- sapply(load.libs,require,character.only = TRUE)
if(all(status)){
  print("SUCCESS: You have successfully installed and loaded all required libraries.")
} else{
  cat("ERROR: One or more libraries failed to install correctly. Check the following list for FALSE cases and try again...\n\n")
  status
}

bd_genes_fe = read.table("/Users/aak/data/AI-CAN-UK/gdd_v3/ms_gene_var_count_sig.txt",head=T,sep="\t")
bd_genes_fdr = read.table("/Users/aak/data/AI-CAN-UK/gdd_v3/ms_gene_var_count_sig_fdr.txt",head=T,sep="\t")
fdr.genes.ez <- clusterProfiler::bitr(bd_genes_fdr$Gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
x.fdr = enrichPathway(gene.fdr.ez,readable = T)


### PathFindR #####

pf.df.fe = bd_genes_fe[,c("Gene","p_Value_fe")]
pf.df.fdr = bd_genes_fdr[,c("Gene","p_Value_fe")]
colnames(pf.df.fe) = c("Gene_symbol","p_val")
colnames(pf.df.fdr) = c("Gene_symbol","p_val")
custom_result = run_pathfindR(pf.df,
                              gene_sets = "Custom",
                              custom_descriptions = custom_descriptions,
                              max_gset_size = Inf, # DO NOT LIMIT GENE SET SIZE
                              iterations = 1, # for demo, setting number of iterations to 1
                              output_dir = "data/AI-CAN-UK/gdd_v3/BurdenFDR")
                              
                              )

res.pf.fdr.kg = run_pathfindR(pf.df.fdr,
                       adj_method="fdr",
                       gene_sets = c("KEGG"),
                       output_dir = "~/data/AI-CAN-UK/gdd_v3/BurdenFDR_KEGG")

res.pf.fdr.rc = run_pathfindR(pf.df.fdr,
                              adj_method="fdr",
                              gene_sets = c("Reactome"),
                              output_dir = "~/data/AI-CAN-UK/gdd_v3/BurdenFDR_RCT")
