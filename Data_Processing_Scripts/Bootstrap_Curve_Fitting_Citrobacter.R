rm(list=ls())
library(data.table)

calc_lowess<- function(data){
  #fit LOWESS To data and return predicted data.points
  xs = seq(0,1,length=501)
  m <- loess(data$Dissimilarity ~ data$Overlap,span=0.2,family="symmetric",iterations=5,statistics='none')
  y = predict(m,xs)
  return(y)
}

calc_lm<- function(data,med){
  #fit lm to data (above median overlap) and return predicted datapoints.
  xs = seq(0,1,length=501)
  data = data[Overlap>med]
  m <- lm(Dissimilarity ~ Overlap,data)
  y = xs*m$coefficients[2] + m$coefficients[1]
  return(y)
}

calc_lci <- function(sdat,conf){
  #Determine 95% lower CI for curves
  return(apply(sdat,2,function(x)  quantile(na.omit(x),(1-conf)/2)))
}

calc_uci <- function(sdat,conf){
  #Determine 95% upper CI for curves
  return(apply(sdat,2,function(x)  quantile(na.omit(x),1- (1-conf)/2)))
}

frac_pos_slope <- function(sdat){
  #Determine fraction of lm models with mositive slopes
  
  return(sum(sdat[,ncol(sdat)] > sdat [,ncol(sdat)-1])/nrow(sdat)) 
}


save_summary <-function(olm,lowess,fn){
  #Save olm and lowess data with given filename
  xs = seq(0,1,length=501)
  sdata= data.frame(xs,mean_lm=colMeans(olm,na.rm=TRUE),
                    mean_lowess=colMeans(lowess,na.rm=TRUE),
                    UCI_lm = calc_uci(olm,0.95),
                    UCI_lowess = calc_uci(lowess,0.95),
                    LCI_lm = calc_lci(olm,0.95),
                    LCI_lowess = calc_lci(lowess,0.95),
                    PValue = frac_pos_slope(olm))
  fwrite(sdata,fn)
  
  return()
}
generate_curves <- function(dat,n,fn1,fn2){
  #Curve fitting for citrobacter
  odata = fread('../Data/Emergent_Comunity_Data_Equilibrium.csv')[Carbon_Source=='Glucose' & Transfer==12]
  coms = unique(odata[ESV_ID == 'Citrobacter',]$Sample_ID)
  olm<- c()
  lowess<- c()
  no_olm<- c()
  no_lowess<- c()
  for(j in 1:n){
    print(j)
    tdoc_df = dat[Run==j,]
    t<- tdoc_df[Var1 %in% coms & Var2 %in% coms,]
    tNo <- tdoc_df[!c(Var1 %in% coms | Var2 %in%  coms),]
    med1 = median(t$Overlap)
    med2 = median(tNo$Overlap)
    olm<- rbind(olm,calc_lm(t,med1))
    lowess<- rbind(lowess,calc_lowess(t))
    no_olm<- rbind(no_olm,calc_lm(tNo,med2))
    no_lowess <- rbind(no_lowess,calc_lowess(tNo))
  }
  save_summary(olm,lowess,fn1)
  save_summary(no_olm,no_lowess,fn2)
  return()
}


set.seed(1)
null = 0

doc_df_b=fread('../Data/Dissimilarity_Overlap_Bootstrapped.csv')
doc_df_b = doc_df_b[Dissimilarity>0 & Carbon_Source_1 == 'Glucose'& Carbon_Source_2 == 'Glucose',]
generate_curves(doc_df_b,500,'../Stat_Outputs/Curve_Fitting_Citrobacter.csv','../Stat_Outputs/Curve_Fitting_No_Citrobacter.csv')

set.seed(1)
null = 1
doc_df_b=fread('../Data/Dissimilarity_Overlap_Bootstrapped_Null.csv')
doc_df_b = doc_df_b[Dissimilarity>0 & Carbon_Source_1 == 'Glucose'& Carbon_Source_2 == 'Glucose',]
generate_curves(doc_df_b,500,'../Stat_Outputs/Curve_Fitting_Citrobacter_Null.csv','../Stat_Outputs/Curve_Fitting_No_Citrobacter_Null.csv')

