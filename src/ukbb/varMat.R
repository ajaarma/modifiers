#########################################################
# Description: Get the variant Matrix for a given gene
#
# Author: Ajay A. Kumar (aak@ebi.ac.uk)
#
#
########################################################

require('dplyr')
require('Matrix')
getVarMat = function(df,mapData,id_samples){

    message('Build Variant Matrix')
    df = data.frame(df)
    ukb_samples = data.frame(mapData)[,1]

    #row_names = unique(as.vector(df[,2]))
    row_names = ukb_samples[!ukb_samples %in% id_samples]
    col_names = unique(as.vector(df[,1]))

    mat = mat.or.vec(length(row_names),length(col_names))
    colnames(mat) = col_names
    rownames(mat) = row_names

    for (i in 1:dim(df)[1]){
        row_id = df[i,2]
        col_id = df[i,1]
        r_index = which(rownames(mat)==row_id)
        c_index = which(colnames(mat)==col_id)
        mat[r_index,c_index] = 1
    }
    return(mat)
}

mergeGMF <- function(){
    #Subroutine to get GMF data
    load('/hps/nobackup/dunham/ai-uk-can/analysis/matrices/ukbTab.df')
    
    gmfFile = '/nfs/research/dunham/samples/ukbb/data/gmf_test/samples_gmf_test.txt'
    gmfData = fread(gmfFile,sep='\t',header=T,stringsAsFactors=F,quote="")

    gmf.tdi = gmfData[,c('f.eid','f.189.0.0')]
    ukbMerge = dplyr::left_join(ukbDF,gmf.tdi,by=join_by(f.eid_samples==f.eid))

    ind.tdi = which(colnames(ukbMerge)=='f.189.0.0')
    colnames(ukbMerge)[ind.tdi] = 'f.189.0.0_tdi'

    save(ukbMerge, file='/hps/nobackup/dunham/ai-uk-can/analysis/matrices/ukb.cgt.gmf')

    return(ukbMerge)

}

getGmfMat <- function(ukbDF) {
    #Subroutine to get GMF Reg matrices
    

}

getCgtMat <- function(ukbDF,case.df,run_type,out_res_dir){
    #Subroutine to get CGT var matrices
    
    varDF = cbind(rownames(case.df),as.data.frame(as.matrix(case.df)))
    colnames(varDF)[1] = 'f.eid_samples'


    #ukb.var.mat =  dplyr::left_join(varDF,ukbDF,by = join_by(f.eid_samples == f.eid_samples))
    X_alc_df = getAlcoDF(ukbDF)
    X_age_df = getAge(ukbDF)
    X_sex_df = getGender(ukbDF)
 
    if(run_type %in% c("sds","tmtA","tmtB","fit","dst","pmt","rtt")){
        if(run_type == "sds") { #,"tmtA","tmtB","fit","dst","pmt","rtt")){
               
            Y.sds.vec = func.sds.cgt(ukbDF)
            tmpDF = cbind(Y.sds.vec,X_sex_df,X_age_df,X_alc_df)
            tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)
            sds.cov.vars.df = dplyr::left_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
            cgt.rm.df.raw = rmZeroCols(sds.cov.vars.df)
            save(cgt.rm.df.raw,file='/hps/nobackup/dunham/ai-uk-can/analysis/sds/mat/cgt.sds.raw.df')
            cgt.rm.df = subset(cgt.rm.df.raw, (sds >3 & sds <33))
            save(cgt.rm.df.raw,file='/hps/nobackup/dunham/ai-uk-can/analysis/sds/mat/cgt.sds.df')
            return(cgt.rm.df)

        }else if(run_type == 'tmtA'){
            Y.tmtA.vec = func.tmtA.cgt(ukbDF)
            Y.tmtA.log = log(Y.tmtA.vec)
            tmpDF = cbind(Y.tmtA.log,X_sex_df,X_age_df,X_alc_df)
            tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)
            tmtA.cov.vars.df = dplyr::left_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
            tmtA.cov.vars.nna = subset(tmtA.cov.vars.df,!is.na(tmtA))
            cgt.rm.df = rmZeroCols(tmtA.cov.vars.nna)
            save(cgt.rm.df,file='/hps/nobackup/dunham/ai-uk-can/analysis/sds/mat/cgt.tmtA.df')
            return(cgt.rm.df)
        
        }else if(run_type == 'tmtB'){
            Y.tmtB.vec = func.tmtB.cgt(ukbDF)
            Y.tmtB.log = log(Y.tmtB.vec)
            tmpDF = cbind(Y.tmtB.log,X_sex_df,X_age_df,X_alc_df)
            tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)
            tmtB.cov.vars.df = dplyr::left_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
            tmtB.cov.vars.nna = subset(tmtB.cov.vars.df,!is.na(tmtB))
            cgt.rm.df = rmZeroCols(tmtA.cov.vars.nna)
            return(cgt.rm.df)
        
        }else if(run_type == 'fit') {
            Y.fit.vec = func.fit.cgt(ukbDF)
            tmpDF = cbind(Y.fit.vec,X_sex_df,X_age_df,X_alc_df)
            tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)
            fit.cov.vars.df = dplyr::left_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
            fit.cov.vars.nna = subset(fit.cov.vars.df,!is.na(fit))
            cgt.rm.df = rmZeroCols(fit.cov.vars.nna)
            save(cgt.rm.df, file='/hps/nobackup/dunham/ai-uk-can/analysis/fit/mat/cgt.fit.df')
            return(cgt.rm.df)
        
        }else if (run_type == 'dst'){
            Y.dst.vec = func.dst.cgt(ukbDF) 
            tmpDF = cbind(Y.dst.vec,X_sex_df,X_age_df,X_alc_df)
            tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)
            dst.cov.vars.df = dplyr::left_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
            dst.cov.vars.nna = subset(dst.cov.vars.df,!is.na(dst))
            cgt.rm.df = rmZeroCols(dst.cov.vars.nna)
            save(cgt.rm.df, file='/hps/nobackup/dunham/ai-uk-can/analysis/dst/mat/cgt.dst.df')
            return(cgt.rm.df)
        
        }else if (run_type == 'pmt'){
            Y.pmt.vec = func.pmt.cgt(ukbDF) 
            Y.pmt.log = log(1+Y.pmt.vec)
            tmpDF = cbind(Y.pmt.log,X_sex_df,X_age_df,X_alc_df) 
            tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)
            pmt.cov.vars.df = dplyr::left_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
            pmt.cov.vars.nna = subset(pmt.cov.vars.df,!is.na(pmt))
            pmt.cov.vars.nna.nz = subset(pmt.cov.vars.nna,pmt !=0)
            cgt.rm.df = rmZeroCols(pmt.cov.vars.nna.nz)
            save(cgt.rm.df,file='/hps/nobackup/dunham/ai-uk-can/analysis/pmt/mat/cgt.pmt.df')
            
            return(cgt.rm.df)
        }else if (run_type == 'rtt'){
            Y.rt.vec = func.rtt.cgt(ukbDF)
            tmpDF = cbind(Y.rt.vec,X_sex_df,X_age_df,X_alc_df)
            tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)
            rtt.cov.vars.df.1 = dplyr::left_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
            rtt.cov.vars.df.2 = subset(rtt.cov.vars.df.1,rtt>200 & rtt <1200)
            cgt.rm.df = rtt.cov.vars.df.2
            save(cgt.rm.df,file='/hps/nobackup/dunham/ai-uk-can/analysis/pmt/mat/cgt.rtt.df')

        } 
        
    }else if (run_type == 'tdi'){

        tmpDF = cbind(cbind(X_sex_df,X_age_df))
        tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)
        varMat.cov = dplyr::left_join(varDF,tmpDF,by = join_by(f.eid_samples == f.eid_samples))

        message('Extract TDI and covariates related info')
        Y_tdi_vec = getTdiDF(ukbDF,'tdi')
        X_ncars_vec = getTdiDF(ukbDF,'ncars')
        X_nhh_vec = getTdiDF(ukbDF,'nih')
        X_rent_vec = getTdiDF(ukbDF,'rent')
        X_emp_vec = getTdiDF(ukbDF,'emp')

        tdi.comb.df = as.data.frame(cbind(Y_tdi_vec,X_ncars_vec,X_nhh_vec,X_rent_vec,X_emp_vec))
        tdi.comb.df.2 = tdi.comb.df[!is.na(tdi.comb.df$tdi),]

        tdi.ncars.outList = func.tdi.ncars.cor(tdi.comb.df.2)
        tdi.ncars.df = tdi.ncars.outList[[1]]

        tdi.nhh.outList = func.tdi.nhh.cor(tdi.comb.df.2)
        tdi.nhh.df = tdi.nhh.outList[[1]]

        tdi.rent.outList = func.tdi.rent.cor(tdi.comb.df.2)
        tdi.rent.df = tdi.rent.outList[[1]]

        tdi.emp.outList = func.tdi.emp.cor(tdi.comb.df.2)
        tdi.emp.df = tdi.emp.outList[[1]]

        message('Intersection between tdi & covariates')
        tdi.cov.common = Reduce(intersect,list(rownames(Y_tdi_vec), rownames(tdi.ncars.df),
                                               rownames(tdi.nhh.df), rownames(tdi.rent.df),
                                               rownames(tdi.emp.df)
                                               ))
		tdi.vec.df2 = Y_tdi_vec[match(tdi.cov.common,rownames(Y_tdi_vec)),]
		tdi.ncars.df2 = tdi.ncars.df[match(tdi.cov.common,rownames(tdi.ncars.df)),]
		tdi.nhh.df2 = tdi.nhh.df[match(tdi.cov.common,rownames(tdi.nhh.df)),]
		tdi.rent.df2 = tdi.rent.df[match(tdi.cov.common,rownames(tdi.rent.df)),]
		tdi.emp.df2 = tdi.emp.df[match(tdi.cov.common,rownames(tdi.emp.df)),]
		tdi.varMat.cov.df2 = varMat.cov[match(tdi.cov.common,varMat.cov$f.eid_samples),]

		tdi.cov.df = data.frame(cbind(tdi.vec.df2,tdi.ncars.df2$ncars,tdi.nhh.df2$nhh,tdi.rent.df2$rent,tdi.emp.df2$emp))
		colnames(tdi.cov.df) = c("tdi","ncars","nhh","rent","emp")
        rownames(tdi.cov.df) = tdi.cov.common

        tdi.cov.vars.df = cbind(tdi.cov.df,tdi.varMat.cov.df2[,-1])
        
        reg.rm.df = tdi.cov.vars.df

        return(reg.rm.df)
    
    }

}

mergeFrameMissVars <- function(anal_type){

    if (anal_type == 'both') {
        load('/hps/nobackup/dunham/ai-uk-can/analysis/impact_ukbb/frame/frame_impact_list')
        frame.df = outList[[1]]
        load('/hps/nobackup/dunham/ai-uk-can/analysis/impact_ukbb/miss/miss_impact_list')
        miss.df = outList[[1]]

        var.frame.df = cbind(rownames(frame.df),as.data.frame(as.matrix(frame.df)))
        var.miss.df = cbind(rownames(miss.df),as.data.frame(as.matrix(miss.df)))

        colnames(var.miss.df)[1] = 'f.eid_samples'
        colnames(var.frame.df)[1] = 'f.eid_samples'

        fm.merge = dplyr::full_join(var.frame.df,var.miss.df,by='f.eid_samples')
        gene.sum = data.frame(apply(fm.merge[,-1],1,sum)); colnames(gene.sum)[1] = 'gene.sum'
        gene.sum.df = cbind(fm.merge$f.eid_samples,gene.sum);colnames(gene.sum.df)[1] = 'f.eid_samples'
    }else if (anal_type == 'miss'){
        load('/hps/nobackup/dunham/ai-uk-can/analysis/impact_ukbb/miss/miss_impact_list')
        miss.df.raw = outList[[1]]
        miss.df = rmGenesLT10(miss.df.raw,colnames(miss.df.raw))

        var.miss.df = cbind(rownames(miss.df),as.data.frame(as.matrix(miss.df)))
        colnames(var.miss.df)[1] = 'f.eid_samples'
        fm.merge = var.miss.df
        gene.sum = data.frame(apply(fm.merge[,-1],1,sum)); colnames(gene.sum)[1] = 'gene.sum'
        gene.sum.df = cbind(fm.merge$f.eid_samples,gene.sum);colnames(gene.sum.df)[1] = 'f.eid_samples'
    
    }else if (anal_type == 'frame'){
        load('/hps/nobackup/dunham/ai-uk-can/analysis/impact_ukbb/frame/frame_impact_list')
        frame.df.raw = outList[[1]]
        frame.df = rmGenesLT10(frame.df.raw,colnames(frame.df.raw))
        
        var.frame.df = cbind(rownames(frame.df),as.data.frame(as.matrix(frame.df)))
        colnames(var.frame.df)[1] = 'f.eid_samples'
        fm.merge = var.frame.df
        gene.sum = data.frame(apply(fm.merge[,-1],1,sum)); colnames(gene.sum)[1] = 'gene.sum'
        gene.sum.df = cbind(fm.merge$f.eid_samples,gene.sum);colnames(gene.sum.df)[1] = 'f.eid_samples'
    }
    return(gene.sum.df)
}


mergeCov <- function(gene.sum.df,anal_type) {

    message('--Adding other covariates')
    wdSamples = read.table('/nfs/research/dunham/samples/ukbb/data/450k/samplesINpheno_ToBeFiltered.txt',header=F)
    gene.sum.filt = gene.sum.df[!gene.sum.df$f.eid_samples %in% wdSamples$V1,]

    #---------------------------------------------
    # Load UKB - cgt-gmf all white british dataset
    #---------------------------------------------
    load('/hps/nobackup/dunham/ai-uk-can/analysis/matrices/ukb.cgt.gmf.white_brit')
    X_sex_df = getGender(ukbDF)
    X_age_df = getStdAge(ukbDF)
    X_age_sq = X_age_df^2
    colnames(X_age_sq) = "age_sq"
    X_sex_age = X_sex_df[,-1] * X_age_df
    colnames(X_sex_age) = "sex_age"
    X_sex_age_sq = X_sex_df[,-1] * X_age_sq
    colnames(X_sex_age_sq) = "sex_age_sq"
    X_alc_df = getAlcoDF(ukbDF)
    X_pcs_df = getUKBpcs(ukbDF)

    if(anal_type == 'fit'){
        Y.fit.vec = func.fit.cgt(ukbDF)
        tmpDF = cbind(Y.fit.vec,X_sex_df,X_age_df,X_alc_df,X_age_sq,X_sex_age,X_sex_age_sq,X_pcs_df)
        tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)
        fit.cov.vars.df = dplyr::inner_join(gene.sum.df,tmpDF,by = join_by(f.eid_samples == f.eid_samples))
        
        #fit.cov.vars.df$fit = scale(fit.cov.vars.df$fit)
        #fit.cov.vars.df$fit = fit.cov.vars.df$fit

        return(fit.cov.vars.df)
    }else if (anal_type =='nedu'){
        Y.qual.vec = getQualDF(ukbDF)
        Y.qual.vec.n1 = subset(Y.qual.vec, qualf !=-3)
        Y.qual.edu = mapEduYears(Y.qual.vec.n1)
        #Y.qual.edu = subset(Y.qual.vec.n1,nedu >4)
        Y.qual.edu$f.eid_samples = as.character(Y.qual.edu$f.eid_samples)
        tmpDF = cbind(X_sex_df,X_age_df,X_alc_df,X_age_sq,X_sex_age,X_sex_age_sq,X_pcs_df)
        tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)

        tmpDF.qual = dplyr::inner_join(tmpDF,Y.qual.edu,by=join_by(f.eid_samples == f.eid_samples))
        tmpDF.qual$f.eid_samples = as.character(tmpDF.qual$f.eid_samples)

        edu.cov.vars.df.1 = dplyr::inner_join(gene.sum.df,tmpDF.qual,by = join_by(f.eid_samples == f.eid_samples))
        edu.cov.vars.df = edu.cov.vars.df.1[,!colnames(edu.cov.vars.df.1) %in% c('qualf')]
        #edu.cov.vars.df = edu.cov.vars.df.1[,!colnames(edu.cov.vars.df.1) %in% c('nedu')]; colnames(edu.cov.vars.df)[30] = 'nedu'

        return(edu.cov.vars.df)
    } else if (anal_type == 'tdi'){
        tmpDF = cbind(X_sex_df,X_age_df,X_alc_df,X_age_sq,X_sex_age,X_sex_age_sq,X_pcs_df)
        tmpDF$f.eid_samples = as.character(tmpDF$f.eid_samples)
		varMat.cov = dplyr::inner_join(gene.sum.df,tmpDF,by = join_by(f.eid_samples == f.eid_samples))        

        message('Extract TDI and covariates related info')
        Y_tdi_vec = getTdiDF(ukbDF,'tdi')
        X_ncars_vec = getTdiDF(ukbDF,'ncars')
        X_nhh_vec = getTdiDF(ukbDF,'nhh')
        X_rent_vec = getTdiDF(ukbDF,'rent')
        X_emp_vec = getTdiDF(ukbDF,'emp')

        f.eid_samples = as.character(ukbDF$f.eid_samples)
        tdi.comb.df = as.data.frame(cbind(f.eid_samples,Y_tdi_vec,X_ncars_vec,X_nhh_vec,X_rent_vec,X_emp_vec))
        tdi.varMat.cov.1 = dplyr::inner_join(tdi.comb.df,varMat.cov,by = join_by(f.eid_samples == f.eid_samples))
            
        tdi.varMat.cov = tdi.varMat.cov.1[,c(1:6)]
        tdi.comb.df.2 = subset(tdi.varMat.cov,!is.na(tdi))

        tdi.ncars.outList = func.tdi.ncars.cor(tdi.comb.df.2)
        tdi.ncars.df = tdi.ncars.outList[[1]]

        tdi.nhh.outList.cor = func.tdi.nhh.cor(tdi.comb.df.2)[[2]]
        tdi.nhh.outList = func.tdi.nhh.cor(tdi.ncars.df)
        tdi.nhh.df = tdi.nhh.outList[[1]]
        tdi.nhh.df.cor = tdi.nhh.outList[[2]]; rm(tdi.nhh.outList)

        tdi.rent.outList.cor = func.tdi.rent.cor(tdi.comb.df.2)[[2]]
        tdi.rent.outList = func.tdi.rent.cor(tdi.nhh.df)
        tdi.rent.df = tdi.rent.outList[[1]]; rm(tdi.rent.outList)

        tdi.emp.outList.cor = func.tdi.emp.cor(tdi.comb.df.2)[[2]]
        tdi.emp.outList = func.tdi.emp.cor(tdi.rent.df)
        tdi.emp.df = tdi.emp.outList[[1]]; rm(tdi.emp.outList)

        tdi.varMat.cov = dplyr::left_join(tdi.emp.df,varMat.cov,by = join_by(f.eid_samples == f.eid_samples))

        return(tdi.varMat.cov)	

    }
}
