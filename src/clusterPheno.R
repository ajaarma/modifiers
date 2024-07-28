.libPaths('/hps/software/users/dunham/R_lib/4.0.3/')

require('data.table')

pat.gene.df = fread('/nfs/research/dunham/samples/ddd/analysis/mat/table_ID_severity_gdd_gene_ID.txt',
                    sep='\t',head=T,stringsAsFactors=TRUE,quote=F)

pat.hpo.df = fread('/nfs/research/dunham/samples/ddd/analysis/mat/table_ID_severity_hpo_gdd_only_num.txt',
                    sep=' ',head=T,stringsAsFactors=TRUE,quote=F)



pat.gene.mat = as.matrix(data.frame(pat.gene.df)[,2:dim(pat.gene.df)[2]])
rownames(pat.gene.mat) = as.vector(pat.gene.df$FID)
pat.hpo.mat = as.matrix(data.frame(pat.hpo.df)[,2:dim(pat.hpo.df)[2]])
rownames(pat.hpo.mat) = as.vector(pat.hpo.df$FID)

####### Prob. Matrix for P(G|I) & P(I|G) ########
prob.gene_pat.mat = apply(pat.gene.mat,2,function(x) {
                                            y = x/sum(x)
                                            return(y)
                    })
prob.gene_pat.mat[is.nan(prob.gene_pat.mat)] <- 0

prob.pat_gene.mat = t(apply(pat.gene.mat,1,function(x) {
                                            y = x/sum(x)
                                            return(y)
                    }))

prob.pat_gene.mat[is.nan(prob.pat_gene.mat)] <- 0


###### Prob. Matrix for P(H|I) & P(I|H) #######

prob.hpo_pat.mat = apply(pat.hpo.mat,2,function(x) {
                                            y = x/sum(x)
                                            return(y)
                    })

prob.hpo_pat.mat[is.nan(prob.pat_gene.mat)] <- 0

prob.pat_hpo.mat = t(apply(pat.hpo.mat,1,function(x) {
                                            y = x/sum(x)
                                            return(y)
                    }))

prob.pat_hpo.mat[is.nan(prob.pat_hpo.mat)] <- 0

### Prob. of choosing an individual with mild, Mod, Severe, General Pheno ###
### To Do : Get the mild,mod,sev,general mapping from Tanai to add these prob
### vector 

p.mild = 639/dim(pat.gene.df)[1]
p.mod = 1417/dim(pat.gene.df)[1]
p.sev = 1040/dim(pat.gene.df)[1]
p.gen = 3502/dim(pat.gene.df)[1]

prior.prob.cat = as.matrix(mat.or.vec(dim(pat.gene.df)[1],1))
rownames(prior.prob.cat) = rownames(prob.pat_hpo.mat)

cat_list = c("Mild","Moderate","Severe","General")

for (ele in cat_list){
    ind_ele = which(grepl(ele,rownames(prior.prob.cat))==TRUE)
    if (ele=="Mild"){
        p_sub = rep(p.mild,length(ind_ele))
    }else if(ele=="Moderate"){
        p_sub = rep(p.mod,length(ind_ele))
    }else if(ele == "Severe"){
        p_sub = rep(p.sev,length(ind_ele))
    }else if(ele=="General"){
        p_sub = rep(p.gen,length(ind_ele))
    }
    
    prior.prob.cat[ind_ele] = p_sub
}

####### Computing P(G|H) distributin matrix ##########
prob.hpo.pat.prod.mat = apply(prob.hpo_pat.mat,2,function(x){
                                                y = x*as.vector(prior.prob.cat)
                                                return(y)
                                              }
                             )

message('-- Computing P(H|I)*P(I) matrix')
prob.gene.hpo.dot.prod.mat = t(prob.gene_pat.mat)%*% prob.hpo.pat.prod.mat #%*% t(prior.prob.cat) 

message('-- Computing P(G|H) = P(G|I)%*% [P(H|)*P(I)] matrix')
prob.gene_hpo.mat = prob.gene.hpo.dot.prod.mat

message("-- Save the output file")
save(prob.gene_hpo.mat,file="/nfs/research/dunham/samples/ddd/data/gene_list/prob.gene_hpo.mat")
