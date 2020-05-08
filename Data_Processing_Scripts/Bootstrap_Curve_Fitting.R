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
generate_curves <- function(dat,n,null){
  # FITS Lowess and OLM for every bootstrap realization for all of the different subests plotted in the paper (each chunk
  # IS fitted with respect to a one median value)
  #Same Environment broken down by Csources pair  and by same inouculum for figure 1 and S1
  SameEnv_olm = c()
  SameEnv_lowess = c()
  SameEnvSameInoc_olm = c()
  SameEnvSameInoc_lowess= c()
  SameEnvDiffInoc_olm = c()
  SameEnvDiffInoc_lowess= c()
  Glc_olm <- c()
  Glc_lowess<- c()
  Cit_olm<- c()
  Cit_lowess<- c()
  Leu_olm<- c()
  Leu_lowess<- c()

  
  #Different Environment broken down by Csources pair for figure s2 and s3
  DiffEnv_olm = c()
  DiffEnv_lowess = c()
  Glc_Cit_olm <- c()
  Glc_Cit_lowess<- c()
  Glc_Leu_olm<- c()
  Glc_Leu_lowess<- c()
  Cit_Leu_olm<- c()
  Cit_Leu_lowess<- c()
  DiffEnvSameInoc_olm = c()
  DiffEnvSameInoc_lowess = c()
  DiffEnvDiffInoc_olm = c()
  DiffEnvDiffInoc_lowess = c()
  
  #Break down by different combinations of same environment and same inoculum for figure s10
  All_olm<- c() #All
  All_lowess<- c() #All
  All_SameEnvSameInoc_olm = c() 
  All_SameEnvSameInoc_lowess = c()
  All_DiffEnvSameInoc_olm = c()
  All_DiffEnvSameInoc_lowess = c()
  All_SameEnvDiffInoc_olm = c()
  All_SameEnvDiffInoc_lowess = c()
  All_DiffEnvDiffInoc_olm = c()
  All_DiffEnvDiffInoc_lowess = c()

  for(j in 1:n){
    print(j)
    #Data for all pairs
    All = dat[Run==j,]
    
    #First data for same environment
    SameEnvdf = All[Same_Environment == TRUE]
    Glcdf <- SameEnvdf[Carbon_Source_1 == 'Glucose',]
    Citdf <-SameEnvdf[Carbon_Source_1 == 'Citrate',]
    Leudf <- SameEnvdf[Carbon_Source_1 == 'Leucine',]
    SameEnvSameInocdf = SameEnvdf[Same_Inoculum==TRUE]
    SameEnvDiffInocdf = SameEnvdf[Same_Inoculum==FALSE]
    
    #next data for different environment
    DiffEnvdf = All[Same_Environment == FALSE]
    GlcCitdf <- DiffEnvdf[Carbon_Source_1 == 'Glucose' & Carbon_Source_2 == 'Citrate',]
    GlcLeudf <-DiffEnvdf[Carbon_Source_1 == 'Leucine' & Carbon_Source_2 == 'Glucose',]
    CitLeudf <- DiffEnvdf[Carbon_Source_1 == 'Leucine' & Carbon_Source_2 == 'Citrate',]
    DiffEnvSameInocdf = DiffEnvdf[Same_Inoculum == TRUE]
    DiffEnvDiffInocdf = DiffEnvdf[Same_Inoculum == FALSE]
    
    #first fit data for same environment (using median for same environment pairs)
    SameEnv_olm<- rbind(SameEnv_olm,calc_lm(SameEnvdf,median(SameEnvdf$Overlap)))
    SameEnv_lowess<- rbind(SameEnv_lowess,calc_lowess(SameEnvdf))
    SameEnvSameInoc_olm  = rbind(SameEnvSameInoc_olm,calc_lm(SameEnvSameInocdf,median(SameEnvdf$Overlap)))
    SameEnvSameInoc_lowess = rbind(SameEnvSameInoc_lowess,calc_lowess(SameEnvSameInocdf))
    SameEnvDiffInoc_olm = rbind(SameEnvDiffInoc_olm,calc_lm(SameEnvDiffInocdf,median(SameEnvdf$Overlap)))
    SameEnvDiffInoc_lowess = rbind(SameEnvDiffInoc_lowess,calc_lowess(SameEnvDiffInocdf))
    Glc_olm<- rbind(Glc_olm,calc_lm(Glcdf,median(SameEnvdf$Overlap)))
    Glc_lowess<- rbind(Glc_lowess,calc_lowess(Glcdf))
    Cit_olm<- rbind(Cit_olm,calc_lm(Citdf,median(SameEnvdf$Overlap)))
    Cit_lowess<- rbind(Cit_lowess,calc_lowess(Citdf))
    Leu_olm<- rbind(Leu_olm,calc_lm(Leudf,median(SameEnvdf$Overlap)))
    Leu_lowess<- rbind(Leu_lowess,calc_lowess(Leudf))
    
    #next fit data for different environment (using median for different environment pairs)
    DiffEnv_olm<- rbind(DiffEnv_olm,calc_lm(DiffEnvdf,median(DiffEnvdf$Overlap)))
    DiffEnv_lowess<- rbind(DiffEnv_lowess,calc_lowess(DiffEnvdf))
    DiffEnvSameInoc_olm  = rbind(DiffEnvSameInoc_olm,calc_lm(DiffEnvSameInocdf,median(DiffEnvdf$Overlap)))
    DiffEnvSameInoc_lowess = rbind(DiffEnvSameInoc_lowess,calc_lowess(DiffEnvSameInocdf))
    DiffEnvDiffInoc_olm = rbind(DiffEnvDiffInoc_olm,calc_lm(DiffEnvDiffInocdf,median(DiffEnvdf$Overlap)))
    DiffEnvDiffInoc_lowess = rbind(DiffEnvDiffInoc_lowess,calc_lowess(DiffEnvDiffInocdf))
    Glc_Cit_olm<- rbind(Glc_Cit_olm,calc_lm(GlcCitdf,median(DiffEnvdf$Overlap)))
    Glc_Cit_lowess<- rbind(Glc_Cit_lowess,calc_lowess(GlcCitdf))
    Glc_Leu_olm<- rbind(Glc_Leu_olm,calc_lm(GlcLeudf,median(DiffEnvdf$Overlap)))
    Glc_Leu_lowess<- rbind(Glc_Leu_lowess,calc_lowess(GlcLeudf))
    Cit_Leu_olm<- rbind(Cit_Leu_olm,calc_lm(CitLeudf,median(DiffEnvdf$Overlap)))
    Cit_Leu_lowess<- rbind(Cit_Leu_lowess,calc_lowess(CitLeudf))
    
    #Finally different combinations of same environment and same inoculum using all is reference point for figure s10
    All_olm<- rbind(All_olm,calc_lm(All,median(All$Overlap))) #All
    All_lowess<-  rbind(All_lowess,calc_lowess(All)) #All
    All_SameEnvSameInoc_olm =  rbind(All_SameEnvSameInoc_olm,calc_lm(SameEnvSameInocdf,median(All$Overlap)))
    All_SameEnvSameInoc_lowess = rbind(All_SameEnvSameInoc_lowess,calc_lowess(SameEnvSameInocdf))
    All_DiffEnvSameInoc_olm = rbind(All_DiffEnvSameInoc_olm,calc_lm(DiffEnvSameInocdf,median(All$Overlap)))
    All_DiffEnvSameInoc_lowess =  rbind(All_DiffEnvSameInoc_lowess,calc_lowess(DiffEnvSameInocdf))
    All_SameEnvDiffInoc_olm =  rbind(All_SameEnvDiffInoc_olm,calc_lm(SameEnvDiffInocdf,median(All$Overlap)))
    All_SameEnvDiffInoc_lowess =  rbind(All_SameEnvDiffInoc_lowess,calc_lowess(SameEnvDiffInocdf))
    All_DiffEnvDiffInoc_olm = rbind(All_DiffEnvDiffInoc_olm,calc_lm(DiffEnvDiffInocdf,median(All$Overlap)))
    All_DiffEnvDiffInoc_lowess =  rbind(All_DiffEnvDiffInoc_lowess,calc_lowess(DiffEnvDiffInocdf))
  }
  if(null==0){
    #ALL
    save_summary(All_olm,All_lowess,'../Stat_Outputs/Curve_Fitting_All.csv')
    save_summary(All_SameEnvSameInoc_olm,All_SameEnvSameInoc_lowess,'../Stat_Outputs/Curve_Fitting_All_SameEnvSameInoc.csv')
    save_summary(All_DiffEnvSameInoc_olm,All_DiffEnvSameInoc_lowess,'../Stat_Outputs/Curve_Fitting_All_DiffEnvSameInoc.csv')
    save_summary(All_SameEnvDiffInoc_olm,All_SameEnvDiffInoc_lowess,'../Stat_Outputs/Curve_Fitting_All_SameEnvDiffInoc.csv')
    save_summary(All_DiffEnvDiffInoc_olm,All_DiffEnvDiffInoc_lowess,'../Stat_Outputs/Curve_Fitting_All_DiffEnvDiffInoc.csv')

    #SAME ENV
    save_summary(SameEnv_olm,SameEnv_lowess,'../Stat_Outputs/Curve_Fitting_SameEnv.csv')
    save_summary(SameEnvSameInoc_olm,SameEnvSameInoc_lowess,'../Stat_Outputs/Curve_Fitting_SameEnvSameInoc.csv')
    save_summary(SameEnvDiffInoc_olm,SameEnvDiffInoc_lowess,'../Stat_Outputs/Curve_Fitting_SameEnvDiffInoc.csv')
    save_summary(Glc_olm,Glc_lowess,'../Stat_Outputs/Curve_Fitting_Glucose.csv')
    save_summary(Cit_olm,Cit_lowess,'../Stat_Outputs/Curve_Fitting_Citrate.csv')
    save_summary(Leu_olm,Leu_lowess,'../Stat_Outputs/Curve_Fitting_Leucine.csv')
    
    #DIfferent Env
    save_summary(DiffEnv_olm,DiffEnv_lowess,'../Stat_Outputs/Curve_Fitting_DiffEnv.csv')
    save_summary(DiffEnvSameInoc_olm,DiffEnvSameInoc_lowess,'../Stat_Outputs/Curve_Fitting_DiffEnvSameInoc.csv')
    save_summary(DiffEnvDiffInoc_olm,DiffEnvDiffInoc_lowess,'../Stat_Outputs/Curve_Fitting_DiffEnvDiffInoc.csv')
    save_summary(Glc_Cit_olm,Glc_Cit_lowess,'../Stat_Outputs/Curve_Fitting_Glucose_Citrate.csv')
    save_summary(Glc_Leu_olm,Glc_Leu_lowess,'../Stat_Outputs/Curve_Fitting_Glucose_Leucine.csv')
    save_summary(Cit_Leu_olm,Cit_Leu_lowess,'../Stat_Outputs/Curve_Fitting_Citrate_Leucine.csv')
    
    
  } else{
    #ALL
    save_summary(All_olm,All_lowess,'../Stat_Outputs/Curve_Fitting_All_Null.csv')
    save_summary(All_SameEnvSameInoc_olm,All_SameEnvSameInoc_lowess,'../Stat_Outputs/Curve_Fitting_All_SameEnvSameInoc_Null.csv')
    save_summary(All_DiffEnvSameInoc_olm,All_DiffEnvSameInoc_lowess,'../Stat_Outputs/Curve_Fitting_All_DiffEnvSameInoc_Null.csv')
    save_summary(All_SameEnvDiffInoc_olm,All_SameEnvDiffInoc_lowess,'../Stat_Outputs/Curve_Fitting_All_SameEnvDiffInoc_Null.csv')
    save_summary(All_DiffEnvDiffInoc_olm,All_DiffEnvDiffInoc_lowess,'../Stat_Outputs/Curve_Fitting_All_DiffEnvDiffInoc_Null.csv')
    
    #SAME ENV
    save_summary(SameEnv_olm,SameEnv_lowess,'../Stat_Outputs/Curve_Fitting_SameEnv_Null.csv')
    save_summary(SameEnvSameInoc_olm,SameEnvSameInoc_lowess,'../Stat_Outputs/Curve_Fitting_SameEnvSameInoc_Null.csv')
    save_summary(SameEnvDiffInoc_olm,SameEnvDiffInoc_lowess,'../Stat_Outputs/Curve_Fitting_SameEnvDiffInoc_Null.csv')
    save_summary(Glc_olm,Glc_lowess,'../Stat_Outputs/Curve_Fitting_Glucose_Null.csv')
    save_summary(Cit_olm,Cit_lowess,'../Stat_Outputs/Curve_Fitting_Citrate_Null.csv')
    save_summary(Leu_olm,Leu_lowess,'../Stat_Outputs/Curve_Fitting_Leucine_Null.csv')
    
    #DIfferent Env
    save_summary(DiffEnv_olm,DiffEnv_lowess,'../Stat_Outputs/Curve_Fitting_DiffEnv_Null.csv')
    save_summary(DiffEnvSameInoc_olm,DiffEnvSameInoc_lowess,'../Stat_Outputs/Curve_Fitting_DiffEnvSameInoc_Null.csv')
    save_summary(DiffEnvDiffInoc_olm,DiffEnvDiffInoc_lowess,'../Stat_Outputs/Curve_Fitting_DiffEnvDiffInoc_Null.csv')
    save_summary(Glc_Cit_olm,Glc_Cit_lowess,'../Stat_Outputs/Curve_Fitting_Glucose_Citrate_Null.csv')
    save_summary(Glc_Leu_olm,Glc_Leu_lowess,'../Stat_Outputs/Curve_Fitting_Glucose_Leucine_Null.csv')
    save_summary(Cit_Leu_olm,Cit_Leu_lowess,'../Stat_Outputs/Curve_Fitting_Citrate_Leucine_Null.csv')

  }
  return()
}


set.seed(1)
doc_df_b=fread('../Data/Dissimilarity_Overlap_Bootstrapped.csv')
doc_df_b = doc_df_b[Dissimilarity>0,]
generate_curves(doc_df_b,500,0)
set.seed(1)
doc_df_b=fread('../Data/Dissimilarity_Overlap_Bootstrapped_Null.csv')
doc_df_b = doc_df_b[Dissimilarity>0,]
generate_curves(doc_df_b,500,1)
