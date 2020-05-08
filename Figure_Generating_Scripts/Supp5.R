rm(list=ls())
library(data.table)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(ggpubr)
library(stringr)
library(stats)
library(zoo)
library(grid)
library(lme4)

set.seed(1)

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

plot_comp_diff <- function(doc_df,o,fn){
  Ad <- doc_df[Carbon_Source_1 == 'Glucose' & Same_Environment == TRUE & Overlap>o,]
  Bd <- doc_df[Carbon_Source_1 == 'Citrate' &Same_Environment == TRUE & Overlap>o,]
  Cd <- doc_df[Carbon_Source_2 == 'Citrate'& Carbon_Source_1 =='Glucose' & Overlap>o,]
  A = rbind(Ad,Bd,Cd)
  A$Treatment = 'NA'
  A[Same_Environment == TRUE & Carbon_Source_1 == 'Glucose',]$Treatment = 'A'
  A[Same_Environment == TRUE & Carbon_Source_1 == 'Citrate',]$Treatment = 'B'
  A[Same_Environment == FALSE,]$Treatment ='C'
  
  #Replace P Values by bootstrapped data
  tests = fread(fn)
  tests = tests[Threshold ==o & !is.na(tests$t)]
  pvals = data.frame()
  for(x in unique(tests$Comparison)){
    t = mean(tests[Comparison==x]$t)
    if(t <0){
      pvals = rbind(pvals,data.frame(Comparison = x,p = sum(tests[Comparison == x]$t>0)/nrow(tests[Comparison == x])))}
    else{
      pvals = rbind(pvals,data.frame(Comparison = x,p = sum(tests[Comparison == x]$t<0)/nrow(tests[Comparison == x])))
    }
  }
  n = length(unique(tests$Run))
  pvals$adj = paste('p < ', ceiling(1000/n)/1000)
  pvals[pvals$p !=0,]$adj = as.character(floor(pvals[pvals$p !=0,]$p*1000)/1000)
  pvals[pvals$p!=0]
  pvals$group1 = c('A','B','A')
  pvals$group2 = c('B','C','C')
  
  p1 <- ggboxplot(A,x='Treatment',y='Dissimilarity',col='Treatment',palette = "jco",
                  add = "jitter",legend='right') +
    # stat_compare_means(comparisons = mycomparisons,method='t.test',size=2) + scale_y_continuous(breaks=c(0,0.8)) +
    guides(col=FALSE) +   labs(x = '') +
    stat_pvalue_manual(pvals,label='adj',y.position = c(0.8,0.85,0.9),size=3) +
    theme(legend.position = '' ,
          axis.line = element_line(size=1),
          axis.text = element_text(size=8)) +    
    scale_x_discrete(labels=c(expression('Glc-Glc'),expression('Cit-Cit'),expression('Glc-Cit')))  +
    ggtitle(paste('t =', o))
  return(p1)
}

doc_df = fread('../Data/Dissimilarity_Overlap.csv')
doc_df = doc_df[Dissimilarity!=0,]
p1<- plot_comp_diff(doc_df,0.9,fn= '../Stat_Outputs/TTest_CSource.csv')
p2<- plot_comp_diff(doc_df,0.91,fn= '../Stat_Outputs/TTest_CSource.csv')
p3<- plot_comp_diff(doc_df,0.92,fn= '../Stat_Outputs/TTest_CSource.csv')
p4<- plot_comp_diff(doc_df,0.93,fn= '../Stat_Outputs/TTest_CSource.csv')
p5<- plot_comp_diff(doc_df,0.94,fn= '../Stat_Outputs/TTest_CSource.csv')
p6<- plot_comp_diff(doc_df,0.95,fn= '../Stat_Outputs/TTest_CSource.csv')
p7<- plot_comp_diff(doc_df,0.96,fn= '../Stat_Outputs/TTest_CSource.csv')
p8<- plot_comp_diff(doc_df,0.97,fn= '../Stat_Outputs/TTest_CSource.csv')
p9<- plot_comp_diff(doc_df,0.98,fn= '../Stat_Outputs/TTest_CSource.csv')
p10<- plot_comp_diff(doc_df,0.99,fn= '../Stat_Outputs/TTest_CSource.csv')
p11<- plot_comp_diff(doc_df,0.9925,fn= '../Stat_Outputs/TTest_CSource.csv')
p12<- plot_comp_diff(doc_df,0.995,fn= '../Stat_Outputs/TTest_CSource.csv')

# p = grid.arrange(pA,pB,layout_matrix=rbind(c(1,1,1,1,1,2,2)))

ggsave('../Final_Figures/Supp5.png',grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol=3),height=12,width=9)
