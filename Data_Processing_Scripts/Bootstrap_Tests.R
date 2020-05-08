rm(list=ls())
library(data.table)
library(operators)
t_test_csource <- function(dat,j,o){
  dat = dat[Run ==j & Overlap>o]
  GlcGlc = dat[Carbon_Source_1 =='Glucose' & Carbon_Source_2 == 'Glucose']
  CitCit = dat[Carbon_Source_1 == 'Citrate' & Carbon_Source_2 == 'Citrate']
  GlcCit = dat[Carbon_Source_1 == 'Glucose' & Carbon_Source_2 == 'Citrate']
  t1  =t.test(GlcGlc$Dissimilarity,CitCit$Dissimilarity)
  t2 = t.test(GlcGlc$Dissimilarity,GlcCit$Dissimilarity)
  t3 = t.test(CitCit$Dissimilarity,GlcCit$Dissimilarity)
  
  return(data.frame(Comparison = c('GlcGlc-CitCit','GlcGlc-GlcCit','CitCit-GlcCit'),
    Threshold =o,
    Run =j,
    t = c(t1$statistic,t2$statistic,t3$statistic),
    p = c(t1$p.value,t2$p.value,t3$p.value)))
}

t_test_citrobacter <- function(dat,j,o){
  odata = fread('../Data/Emergent_Comunity_Data_Equilibrium.csv')[Carbon_Source=='Glucose' & Transfer==12]
  cit_coms = unique(odata[ESV_ID == 'Citrobacter',]$Sample_ID)
  dat = dat[Run ==j & Overlap>o]
  dat$Citrobacter <- 'One'
  dat[Var1 %!in% cit_coms & Var2 %!in% cit_coms,]$Citrobacter= 'Neither'
  dat[Var1 %in% cit_coms & Var2 %in% cit_coms,]$Citrobacter= 'Both'
  dat$Citrobacter = factor(dat$Citrobacter,levels=c('Both','One','Neither'))
  Both = dat[Citrobacter =='Both',]
  One = dat[Citrobacter == 'One',]
  Neither = dat[Citrobacter == 'Neither']
  t1  =t.test(Both$Dissimilarity,One$Dissimilarity)
  t2 = t.test(Both$Dissimilarity,Neither$Dissimilarity)
  t3 = t.test(One$Dissimilarity,Neither$Dissimilarity)
  return(data.frame(Comparison = c('Both-One','Both-Neither','One-Neither'),
                    Threshold =o,
                    Run =j,
                    t = c(t1$statistic,t2$statistic,t3$statistic),
                    p = c(t1$p.value,t2$p.value,t3$p.value)))
}

t_test_raoultella <- function(dat,j,o){
  odata = fread('../Data/Emergent_Comunity_Data_Equilibrium.csv')[Carbon_Source=='Citrate' & Transfer==12]
  cit_coms = unique(odata[ESV_ID == 'Raoultella',]$Sample_ID)
  dat = dat[Run ==j & Overlap>o]
  dat$Citrobacter <- 'One'
  dat[Var1 %!in% cit_coms & Var2 %!in% cit_coms,]$Citrobacter= 'Neither'
  dat[Var1 %in% cit_coms & Var2 %in% cit_coms,]$Citrobacter= 'Both'
  dat$Citrobacter = factor(dat$Citrobacter,levels=c('Both','One','Neither'))
  Both = dat[Citrobacter =='Both',]
  One = dat[Citrobacter == 'One',]
  Neither = dat[Citrobacter == 'Neither']
  t1  =t.test(Both$Dissimilarity,One$Dissimilarity)
  t2 = t.test(Both$Dissimilarity,Neither$Dissimilarity)
  t3 = t.test(One$Dissimilarity,Neither$Dissimilarity)
  return(data.frame(Comparison = c('Both-One','Both-Neither','One-Neither'),
                    Threshold =o,
                    Run =j,
                    t = c(t1$statistic,t2$statistic,t3$statistic),
                    p = c(t1$p.value,t2$p.value,t3$p.value)))
}

t_test_pseudomonas <- function(dat,j,o){
  odata = fread('../Data/Emergent_Comunity_Data_Equilibrium.csv')[Carbon_Source=='Leucine' & Transfer==12]
  cit_coms = unique(odata[ESV_ID == 'Pseudomonas',]$Sample_ID)
  dat = dat[Run ==j & Overlap>o]
  dat$Citrobacter <- 'One'
  dat[Var1 %!in% cit_coms & Var2 %!in% cit_coms,]$Citrobacter= 'Neither'
  dat[Var1 %in% cit_coms & Var2 %in% cit_coms,]$Citrobacter= 'Both'
  dat$Citrobacter = factor(dat$Citrobacter,levels=c('Both','One','Neither'))
  Both = dat[Citrobacter =='Both',]
  One = dat[Citrobacter == 'One',]
  Neither = dat[Citrobacter == 'Neither']
  t1  =t.test(Both$Dissimilarity,One$Dissimilarity)
  t2 = t.test(Both$Dissimilarity,Neither$Dissimilarity)
  t3 = t.test(One$Dissimilarity,Neither$Dissimilarity)
  return(data.frame(Comparison = c('Both-One','Both-Neither','One-Neither'),
                    Threshold =o,
                    Run =j,
                    t = c(t1$statistic,t2$statistic,t3$statistic),
                    p = c(t1$p.value,t2$p.value,t3$p.value)))
}

set.seed(1)
doc_df=fread('../Data/Dissimilarity_Overlap.csv')
doc_df_b=fread('../Data/Dissimilarity_Overlap_Bootstrapped.csv')
doc_df_b_null=fread('../Data/Dissimilarity_Overlap_Bootstrapped_Null.csv')
doc_df = doc_df[Dissimilarity>0,]
doc_df_b = doc_df_b[Dissimilarity>0,]
doc_df_b_null = doc_df_b_null[Dissimilarity>0,]

all_df = data.frame()
all_df_null = data.frame()
same_inoc_df = data.frame()
#Above Median Overlap
cit_df = data.frame()
raoul_df = data.frame()
pseud_df = data.frame()

#Full range of overlaps
cit_df2 = data.frame()

for(j in 1:max(doc_df_b$Run)){
  print(j)
  cit_df = rbind(cit_df,
                 tryCatch(t_test_citrobacter(doc_df_b[Same_Inoculum==TRUE & 
                                                        Carbon_Source_1 =='Glucose' & 
                                                        Carbon_Source_2 == 'Glucose'],j,
                                    median(doc_df[Same_Inoculum==TRUE & 
                                                    Carbon_Source_1 =='Glucose' & 
                                                    Carbon_Source_2 == 'Glucose']$Overlap)),
                          error=function(e){return(data.frame(Comparison = c('Both-One','Both-Neither','One-Neither'),
                                                             Threshold =median(doc_df[Same_Inoculum==TRUE & 
                                                                                        Carbon_Source_1 =='Glucose' & 
                                                                                        Carbon_Source_2 == 'Glucose']$Overlap),
                                                             Run =j,
                                                             t = NA,
                                                             p = NA))})) #This is v tedious but sometimes there are no high overlap pairs in one group in a realization so we just set it as
  raoul_df = rbind(raoul_df,
                 tryCatch(t_test_raoultella(doc_df_b[Same_Inoculum==TRUE & 
                                                        Carbon_Source_1 =='Citrate' & 
                                                        Carbon_Source_2 == 'Citrate'],j,
                                             median(doc_df[Same_Inoculum==TRUE & 
                                                             Carbon_Source_1 =='Citrate' & 
                                                             Carbon_Source_2 == 'Citrate']$Overlap)),
                          error=function(e){return(data.frame(Comparison = c('Both-One','Both-Neither','One-Neither'),
                                                              Threshold =median(doc_df[Same_Inoculum==TRUE & 
                                                                                         Carbon_Source_1 =='Citrate' & 
                                                                                         Carbon_Source_2 == 'Citrate']$Overlap),
                                                              Run =j,
                                                              t = NA,
                                                              p = NA))})) #This is v tedious but sometimes there are no high overlap pairs in one group in a realization so we just set it as
  
  pseud_df = rbind(pseud_df,
                 tryCatch(t_test_pseudomonas(doc_df_b[Same_Inoculum==TRUE & 
                                                        Carbon_Source_1 =='Leucine' & 
                                                        Carbon_Source_2 == 'Leucine'],j,
                                             median(doc_df[Same_Inoculum==TRUE & 
                                                             Carbon_Source_1 =='Leucine' & 
                                                             Carbon_Source_2 == 'Leucine']$Overlap)),
                          error=function(e){return(data.frame(Comparison = c('Both-One','Both-Neither','One-Neither'),
                                                              Threshold = median(doc_df[Same_Inoculum==TRUE & 
                                                                                          Carbon_Source_1 =='Leucine' & 
                                                                                          Carbon_Source_2 == 'Leucine']$Overlap),
                                                              Run =j,
                                                              t = NA,
                                                              p = NA))})) #This is v tedious but sometimes there are no high overlap pairs in one group in a realization so we just set it as
  all_df_null = rbind(all_df_null,t_test_csource(doc_df_b_null,j,0.98))
  same_inoc_df = rbind(same_inoc_df,t_test_csource(doc_df_b[Same_Inoculum==TRUE],j,0.98))
  for(o in c(0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.9925,0.995)){
    #For figures S5 and S6
    all_df = rbind(all_df,tryCatch(t_test_csource(doc_df_b,j,o),
                                   error=function(e){return(data.frame(Comparison = c('GlcGlc-CitCit','GlcGlc-GlcCit','CitCit-GlcCit'),
                                                                                                    Threshold =o,
                                                                                                    Run =j,
                                                                                                    t = NA,
                                                                                                    p = NA))}))
   cit_df2 = rbind(cit_df2,
                   tryCatch(t_test_citrobacter(doc_df_b[Same_Inoculum==TRUE & 
                                                        Carbon_Source_1 =='Glucose' & 
                                                        Carbon_Source_2 == 'Glucose'],j,o),
                          error=function(e){return(data.frame(Comparison = c('Both-One','Both-Neither','One-Neither'),
                                                              Threshold =o,
                                                              Run =j,
                                                              t = NA,
                                                              p = NA))}))    
  }
}
fwrite(all_df,'../Stat_Outputs/TTest_CSource.csv')
fwrite(all_df_null,'../Stat_Outputs/TTest_CSource_Null.csv')
fwrite(cit_df,'../Stat_Outputs/TTest_CSource_Cit.csv')
fwrite(raoul_df,'../Stat_Outputs/TTest_CSource_Raoul.csv')
fwrite(pseud_df,'../Stat_Outputs/TTest_CSource_Pseud.csv')
fwrite(cit_df2,'../Stat_Outputs/TTest_CSource_Cit_Overlap.csv')

