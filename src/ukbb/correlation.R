########################################################
#
# Description: Functions to find correlation between TDI 
#              and other measures of functioning.
#
# Author: Ajay A. Kumar (aak@ebi.ac.uk)
#
#
##########################################################

require('ltm')
require('ggplot2')

func.tdi.ncars.cor <- function(tdi.comb.df.2) {

    #Correlation between TDI and Num. of Cars

    message('Having No cars')
    td.none.df = tdi.comb.df.2
    td.none.df.2 = subset(td.none.df, ncars!=-1) #Remove -1: Do not know
    td.none.df.3 = subset(td.none.df.2, ncars!=-3) #Remove -3: Prefer not to answer
    td.none.df.4 = subset(td.none.df.3, ncars!=0) #Remove 0: NA values orignally for which 0 was substituted
    td.none.df.4[!td.none.df.4$ncars==1,'ncars'] = 0 #Add labels 0: For having at least a car
    td.none.df.4[td.none.df.4$ncars==1,'ncars'] = 1 #Add label 1: For having no car

    bis.cor.ncars = biserial.cor(y=td.none.df.4$ncars,x=td.none.df.4$tdi)
    df = data.frame(cbind(td.none.df.4$tdi,td.none.df.4$ncars))
    colnames(df) = c('tdi','ncars')
    levels(df$ncars) = c("With cars","No Cars")
    m = ggplot(df,aes(y=tdi,x=factor(ncars)))+geom_boxplot()

    #plot(density(td.none.df$tdi,bw=1/50000),
    #    main='Distribution of TDI for UKBB-50K',
    #    xlab='TDI range',ylab='Num. Of Individuals')

    out_list = list()
    out_list[[1]] = td.none.df.4
    out_list[[2]] = bis.cor.ncars

    return(out_list)
}


func.tdi.nhh.cor <- function(tdi.comb.df.2){

    message('Num people in HH')
    td.none.df = tdi.comb.df.2
    td.none.df.2 = subset(td.none.df,nhh!=-1)
    td.none.df.3 = subset(td.none.df.2,nhh!=0)
    td.none.df.4 = subset(td.none.df.3,nhh!=-3)
    td.none.df.4[td.none.df.4$nhh<7,'nhh'] = 0 #Add labels 0: For less overcrowded household
    td.none.df.4[td.none.df.4$nhh>=7,'nhh'] = 1 #Add label 1: For overcorwded household

    bis.cor.nhh = biserial.cor(y=td.none.df.4$nhh,x=td.none.df.4$tdi)
    df = data.frame(cbind(td.none.df.4$tdi,td.none.df.4$nhh))
    colnames(df) = c('tdi','nhh')
    levels(df$nhh) = c("Normal","Overcrowded")
    m = ggplot(df,aes(y=tdi,x=factor(nhh)))+geom_boxplot()

    #plot(density(td.none.df$tdi,bw=1/50000),
    #    main='Distribution of TDI for UKBB-50K',
    #    xlab='TDI range',ylab='Num. Of Individuals')
    out_list = list()
    out_list[[1]] = td.none.df.4
    out_list[[2]] = bis.cor.nhh
    
    return(out_list)
}


func.tdi.rent.cor <- function(tdi.comb.df.2) {

    message('Rent')
    td.none.df = tdi.comb.df.2
    td.none.df.2 = subset(td.none.df,rent!=-7)
    td.none.df.3 = subset(td.none.df.2,rent!=0)
    td.none.df.4 = subset(td.none.df.3,rent!=-3)
    td.none.df.4[td.none.df.4$rent==1,'rent'] = 0 #Add labels 0: For rent other than class 3&4
    td.none.df.4[td.none.df.4$rent==2,'rent'] = 0 #Add labels 0: For rent other than class 3&4
    td.none.df.4[td.none.df.4$rent==3,'rent'] = 1 #Add labels 0: For rent other than class 3&4
    td.none.df.4[td.none.df.4$rent==4,'rent'] = 1 #Add label 1: For rent other than class 3 & 4
    td.none.df.4[td.none.df.4$rent==5,'rent'] = 0 #Add labels 0: For rent from local authority
    td.none.df.4[td.none.df.4$rent==6,'rent'] = 0 #Add label 1: For rent from private


    bis.cor.rent = biserial.cor(y=td.none.df.4$rent,x=td.none.df.4$tdi)
    df = data.frame(cbind(td.none.df.4$tdi,td.none.df.4$rent))
    colnames(df) = c('tdi','rent')
    levels(df$rent) = c("Self-Financed-HH","RentedHH")
    m = ggplot(df,aes(y=tdi,x=factor(rent)))+geom_boxplot()

    out_list = list()
    out_list[[1]] = td.none.df.4
    out_list[[2]] = bis.cor.rent

    return(out_list)
}

func.tdi.emp.cor <- function(tdi.comb.df.2) {
    message('Emp Status')
    td.none.df = tdi.comb.df.2
    td.none.df.2 = subset(td.none.df,emp!=-7)
    td.none.df.3 = subset(td.none.df.2,emp!=0)
    td.none.df.4 = subset(td.none.df.3,emp!=-3)
    td.none.df.4[td.none.df.4$emp==1,'emp'] = 0 #Add labels 0: For rent other than class 3&4
    td.none.df.4[td.none.df.4$emp==2,'emp'] = 0 #Add labels 0: For rent other than class 3&4
    td.none.df.4[td.none.df.4$emp==3,'emp'] = 0 #Add labels 0: For rent other than class 3&4
    td.none.df.4[td.none.df.4$emp==4,'emp'] = 0 #Add label 1: For rent other than class 3 & 4
    td.none.df.4[td.none.df.4$emp==5,'emp'] = 1 #Add labels 0: For rent from local authority
    #td.none.df.4[td.none.df.4$emp==6,'emp'] = 1 #Add label 1: For rent from private
    td.none.df.4[td.none.df.4$emp==6,'emp'] = 0 #Add label 1: For rent from private
    td.none.df.4[td.none.df.4$emp==7,'emp'] = 0 #Add label 1: For rent from private

    bis.cor.emp = biserial.cor(y=td.none.df.4$emp,x=td.none.df.4$tdi)
    df = data.frame(cbind(td.none.df.4$tdi,td.none.df.4$emp))
    colnames(df) = c('tdi','emp')
    levels(df$emp) = c("Employed","Unemployed")
    m = ggplot(df,aes(y=tdi,x=factor(emp)))+geom_boxplot()

    out_list = list()
    out_list[[1]] = td.none.df.4
    out_list[[2]] = bis.cor.emp
    
    return(out_list)
}


