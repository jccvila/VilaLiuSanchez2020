rm(list=ls())
library(data.table)

calc_lowess<- function(data){
  xs = seq(0,1,length=501)
  # m <- try(expr = loess(data$Dissimilarity ~ data$Overlap,span=0.2,family="symmetric",iterations=5,statistics='none'),silent=TRUE)
  # if(class(m) == "try-error"){
  m <- loess(data$Dissimilarity ~ data$Overlap,span=0.2,family="symmetric",iterations=5,statistics='none',surface='direct')
  #   
  # }
  y = predict(m,xs)
  return(y)
}

calc_lm<- function(data,med){
  xs = seq(0,1,length=501)
  data = data[Overlap>med]
  m <- lm(Dissimilarity ~ Overlap,data)
  y = xs*m$coefficients[2] + m$coefficients[1]
  return(y)
}

calc_lci <- function(sdat,conf){
  return(apply(sdat,2,function(x)  quantile(na.omit(x),(1-conf)/2)))
}

calc_uci <- function(sdat,conf){
  return(apply(sdat,2,function(x)  quantile(na.omit(x),1- (1-conf)/2)))
}

frac_pos_slope <- function(sdat){
  return(sum(sdat[,ncol(sdat)] > sdat [,ncol(sdat)-1])/nrow(sdat)) 
}

calc_frac <-function(data){
  return(nrow(data[Overlap>0.5 & Dissimilarity<sqrt(log(2))/2])/nrow(data))
}

is_pos_slope =function(data,med){
  y = calc_lm(data,med)
  val = y[length(y)] > y[length(y)-1]
  return(val)
}

save_summary <-function(olm,lowess,fn){
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
generate_curves(doc_df_b,500,'../Stat_Outputs/Cit_Dat.csv','../Stat_Outputs/No_Cit_Dat.csv')

set.seed(1)
null = 1
doc_df_b=fread('../Data/Dissimilarity_Overlap_Bootstrapped_Null.csv')
doc_df_b = doc_df_b[Dissimilarity>0 & Carbon_Source_1 == 'Glucose'& Carbon_Source_2 == 'Glucose',]
generate_curves(doc_df_b,500,'../Stat_Outputs/Cit_Dat_Null.csv','../Stat_Outputs/No_Cit_Dat_Null.csv')

