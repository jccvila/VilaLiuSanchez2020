rm(list=ls())
library(data.table)

calc_lowess<- function(data){
  xs = seq(0,1,length=501)
  m <- loess(data$Dissimilarity ~ data$Overlap,span=0.2,family="symmetric",iterations=5,statistics='none')
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
generate_curves <- function(dat,n,null){
  # FITS Lowess and OLM for every bootstrap realization for Fig 2A-E
  #Same Environment broken down by Csources
  olm<- c()
  lowess<- c()
  A_olm <- c()
  A_lowess<- c()
  B_olm<- c()
  B_lowess<- c()
  C_olm<- c()
  C_lowess<- c()
  
  #Different Environment broken donw by csource
  D_olm<- c()
  D_lowess<- c()
  E_olm<- c()
  E_lowess<- c()
  F_olm<- c()
  F_lowess<- c()
  G_olm<- c()
  G_lowess<- c()
  
  #Same/Different Environment broken down Inoculum
  H_olm<- c()
  H_lowess<- c()
  I_olm<- c()
  I_lowess<- c()
  J_olm<- c()
  J_lowess<- c()
  K_olm<- c()
  K_lowess<- c()
  
  #All Environment broken down same environment and same inouclum
  L_olm<- c()
  L_lowess<- c()
  M_olm<- c()
  M_lowess<- c()
  N_olm<- c()
  N_lowess<- c()
  O_olm<- c()
  O_lowess<- c()
  P_olm<- c()
  P_lowess<- c()

  for(j in 1:n){
    print(j)
    #First data for same environment (using median for same environment)
    tdoc_df = dat[Run==j & Same_Environment == TRUE,]
    tAd <- tdoc_df[Carbon_Source_1 == 'Glucose',]
    tBd <-tdoc_df[Carbon_Source_1 == 'Citrate',]
    tCd <- tdoc_df[Carbon_Source_1 == 'Leucine',]
    olm<- rbind(olm,calc_lm(tdoc_df,median(tdoc_df$Overlap)))
    lowess<- rbind(lowess,calc_lowess(tdoc_df))
    A_olm<- rbind(A_olm,calc_lm(tAd,median(tdoc_df$Overlap)))
    A_lowess<- rbind(A_lowess,calc_lowess(tAd))
    B_olm<- rbind(B_olm,calc_lm(tBd,median(tdoc_df$Overlap)))
    B_lowess<- rbind(B_lowess,calc_lowess(tBd))
    C_olm<- rbind(C_olm,calc_lm(tCd,median(tdoc_df$Overlap)))
    C_lowess<- rbind(C_lowess,calc_lowess(tCd))
    
    #next data for different environment (using median for different environment)
    tDd = dat[Run==j & Same_Environment == FALSE,]
    tEd <- tDd[Carbon_Source_1 == 'Glucose' & Carbon_Source_2 == 'Citrate',]
    tFd <-tDd[Carbon_Source_1 == 'Leucine' & Carbon_Source_2 == 'Citrate',]
    tGd <- tDd[Carbon_Source_1 == 'Leucine' & Carbon_Source_2 == 'Glucose',]
    D_olm<- rbind(D_olm,calc_lm(tDd,median(tDd$Overlap)))
    D_lowess<- rbind(D_lowess,calc_lowess(tDd))
    E_olm<- rbind(E_olm,calc_lm(tEd,median(tDd$Overlap)))
    E_lowess<- rbind(E_lowess,calc_lowess(tEd))
    F_olm<- rbind(F_olm,calc_lm(tFd,median(tDd$Overlap)))
    F_lowess<- rbind(F_lowess,calc_lowess(tFd))
    G_olm<- rbind(G_olm,calc_lm(tGd,median(tDd$Overlap)))
    G_lowess<- rbind(G_lowess,calc_lowess(tGd))
    
    #next data for all environment (using median for all environment pairs)
    tHd = dat[Run==j & Same_Environment == TRUE & Same_Inoculum ==TRUE,]
    tId = dat[Run==j & Same_Environment == TRUE  & Same_Inoculum ==FALSE,]
    tJd = dat[Run==j & Same_Environment == FALSE  & Same_Inoculum ==TRUE,]
    tKd = dat[Run==j & Same_Environment == FALSE &  Same_Inoculum ==FALSE,]
    tLd = dat[Run==j,]
    
    #Variants for figure 1 by same inouclum
    H_olm<- rbind(H_olm,calc_lm(tHd,median(tdoc_df$Overlap)))
    H_lowess<- rbind(H_lowess,calc_lowess(tHd))
    I_olm<- rbind(I_olm,calc_lm(tId,median(tdoc_df$Overlap)))
    I_lowess<- rbind(I_lowess,calc_lowess(tId))
    
    #Variants for figure S2 by same inocuulm
    J_olm<- rbind(J_olm,calc_lm(tJd,median(tDd$Overlap)))
    J_lowess<- rbind(J_lowess,calc_lowess(tJd))
    K_olm<- rbind(K_olm,calc_lm(tKd,median(tDd$Overlap)))
    K_lowess<- rbind(K_lowess,calc_lowess(tKd))
    
    #Variants for figure S10 (reference global median.)
    L_olm<- rbind(L_olm,calc_lm(tLd,median(tLd$Overlap)))
    L_lowess<- rbind(L_lowess,calc_lowess(tLd))   
    M_olm<- rbind(M_olm,calc_lm(tHd,median(tLd$Overlap)))
    M_lowess<- rbind(M_lowess,calc_lowess(tHd))   
    N_olm<- rbind(N_olm,calc_lm(tId,median(tLd$Overlap)))
    N_lowess<- rbind(N_lowess,calc_lowess(tId))   
    O_olm<- rbind(O_olm,calc_lm(tJd,median(tLd$Overlap)))
    O_lowess<- rbind(O_lowess,calc_lowess(tJd))   
    P_olm<- rbind(P_olm,calc_lm(tKd,median(tLd$Overlap)))
    P_lowess<- rbind(P_lowess,calc_lowess(tKd))   
  }
  if(null==0){
    save_summary(olm,lowess,'../Stat_Outputs/DOC_stats.csv')
    save_summary(A_olm,A_lowess,'../Stat_Outputs/DOC_stats_A.csv')
    save_summary(B_olm,B_lowess,'../Stat_Outputs/DOC_stats_B.csv')
    save_summary(C_olm,C_lowess,'../Stat_Outputs/DOC_stats_C.csv')
    save_summary(D_olm,D_lowess,'../Stat_Outputs/DOC_stats_D.csv')
    save_summary(E_olm,E_lowess,'../Stat_Outputs/DOC_stats_E.csv')
    save_summary(F_olm,F_lowess,'../Stat_Outputs/DOC_stats_F.csv')
    save_summary(G_olm,F_lowess,'../Stat_Outputs/DOC_stats_G.csv')
    save_summary(H_olm,H_lowess,'../Stat_Outputs/DOC_stats_H.csv')
    save_summary(I_olm,I_lowess,'../Stat_Outputs/DOC_stats_I.csv')
    save_summary(J_olm,J_lowess,'../Stat_Outputs/DOC_stats_J.csv')
    save_summary(K_olm,K_lowess,'../Stat_Outputs/DOC_stats_K.csv')    
    save_summary(L_olm,L_lowess,'../Stat_Outputs/DOC_stats_L.csv')    
    save_summary(M_olm,M_lowess,'../Stat_Outputs/DOC_stats_M.csv')    
    save_summary(N_olm,N_lowess,'../Stat_Outputs/DOC_stats_N.csv')    
    save_summary(O_olm,O_lowess,'../Stat_Outputs/DOC_stats_O.csv')    
    save_summary(P_olm,P_lowess,'../Stat_Outputs/DOC_stats_P.csv')    

    
  } else{
    save_summary(olm,lowess,'../Stat_Outputs/DOC_stats_Null.csv')
    save_summary(A_olm,A_lowess,'../Stat_Outputs/DOC_stats_A_Null.csv')
    save_summary(B_olm,B_lowess,'../Stat_Outputs/DOC_stats_B_Null.csv')
    save_summary(C_olm,C_lowess,'../Stat_Outputs/DOC_stats_C_Null.csv')
    save_summary(D_olm,D_lowess,'../Stat_Outputs/DOC_stats_D_Null.csv')
    save_summary(E_olm,E_lowess,'../Stat_Outputs/DOC_stats_E_Null.csv')
    save_summary(F_olm,F_lowess,'../Stat_Outputs/DOC_stats_F_Null.csv')
    save_summary(G_olm,G_lowess,'../Stat_Outputs/DOC_stats_G_Null.csv')
    save_summary(H_olm,H_lowess,'../Stat_Outputs/DOC_stats_H_Null.csv')
    save_summary(I_olm,I_lowess,'../Stat_Outputs/DOC_stats_I_Null.csv')
    save_summary(J_olm,J_lowess,'../Stat_Outputs/DOC_stats_J_Null.csv')
    save_summary(L_olm,L_lowess,'../Stat_Outputs/DOC_stats_L_Null.csv')    
    save_summary(M_olm,M_lowess,'../Stat_Outputs/DOC_stats_M_Null.csv')    
    save_summary(N_olm,N_lowess,'../Stat_Outputs/DOC_stats_N_Null.csv')    
    save_summary(O_olm,O_lowess,'../Stat_Outputs/DOC_stats_O_Null.csv')    
    save_summary(P_olm,P_lowess,'../Stat_Outputs/DOC_stats_P_Null.csv')    
  }
  return()
}
# 
# generate_stats <-function(dat,n,null){
#   # Calculate F and m for every bootstrap realization for figures 2F and supp2
#   f_data = matrix(nrow=9,ncol=n)
#   slope_data = matrix(nrow=9,ncol=n)
#   for(j in 1:n ){
#     print(j)
#     tdoc_df = dat[Run==j,]
#     med = median(tdoc_df$Overlap)
#     Ad <- tdoc_df[Same_Environment == TRUE,]
#     Bd <- tdoc_df[Carbon_Source_1 == 'Glucose'& Same_Environment == TRUE,]
#     Cd <- tdoc_df[Carbon_Source_1 == 'Citrate'& Same_Environment == TRUE,]
#     Dd <- tdoc_df[Carbon_Source_1 == 'Leucine'& Same_Environment == TRUE,]
#     Ed <- tdoc_df[Same_Environment == FALSE,]
#     Fd <- tdoc_df[Carbon_Source_1 == 'Glucose'& Carbon_Source_2 == 'Citrate',]
#     Gd <- tdoc_df[Carbon_Source_1 == 'Leucine'& Carbon_Source_2 == 'Citrate',]
#     Hd <- tdoc_df[Carbon_Source_1 == 'Leucine'& Carbon_Source_2 == 'Glucose',]
#     dat_list = list(tdoc_df,Ad,Bd,Cd,Dd,Ed,Fd,Gd,Hd)
#     f_vec = sapply(dat_list,function(x) calc_frac(x))
#     slope_vec = sapply(dat_list,function(x) is_pos_slope(x,med))
#     f_data[,j] = f_vec
#     slope_data[,j] = slope_vec
#   }
#   f_data = data.frame(f_data)
#   slope_data = data.frame(slope_data)
#   if(null ==0){
#     fwrite(f_data,'../Stat_Outputs/Frac_data.csv')
#     fwrite(slope_data,'../Stat_Outputs/Slope_data.csv')
#   } else{ fwrite(f_data,'../Stat_Outputs/Frac_data_Null.csv')
#     fwrite(slope_data,'../Stat_Outputs/Slope_data_Null.csv')}
#   return()
# }

set.seed(1)
null = 0

if(null == 0){
  doc_df_b=fread('../Data/Dissimilarity_Overlap_Bootstrapped.csv')
}else{
  doc_df_b=fread('../Data/Dissimilarity_Overlap_Bootstrapped_Null.csv')
}
doc_df_b = doc_df_b[Dissimilarity>0,]
generate_curves(doc_df_b,500,null)
set.seed(1)
null = 1

if(null == 0){
  doc_df_b=fread('../Data/Dissimilarity_Overlap_Bootstrapped.csv')
}else{
  doc_df_b=fread('../Data/Dissimilarity_Overlap_Bootstrapped_Null.csv')
}
doc_df_b = doc_df_b[Dissimilarity>0,]
generate_curves(doc_df_b,500,null)
