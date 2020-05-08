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

plot_comp <- function(dat,med,fn){
  odata = fread('../Data/Emergent_Comunity_Data_Equilibrium.csv')[Carbon_Source=='Glucose' & Transfer==12]
  cit_coms = unique(odata[ESV_ID == 'Citrobacter',]$Sample_ID)
  dat$Citrobacter <- 'One'
  dat[Var1 %!in% cit_coms & Var2 %!in% cit_coms,]$Citrobacter= 'Neither'
  dat[Var1 %in% cit_coms & Var2 %in% cit_coms,]$Citrobacter= 'Both'
  dat$Citrobacter = factor(dat$Citrobacter,levels=c('Both','One','Neither'))
  dat = dat[Overlap>med,]
  tests = fread(fn)
  tests = tests[Threshold ==med & !is.na(tests$t)]
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
  pvals$group1 = c('Both','Both','One')
  pvals$group2 = c('One','Neither','Neither')
  pA <- ggboxplot(dat,x='Citrobacter',y='Dissimilarity',col='Citrobacter',palette = "jco",
                  add = "jitter",legend='right') +
    # stat_compare_means(comparisons = mycomparisons,method='t.test',size=2) + scale_y_continuous(breaks=c(0,0.8)) +
    guides(col=FALSE) +   labs(x = '') +
    stat_pvalue_manual(pvals,label='adj',y.position = c(0.8,0.9,0.85),size=3) +
    theme(legend.position = '' ,
          axis.line = element_line(size=1),
          axis.text = element_text(size=8)) +
    ggtitle(paste('t =', med))
  
  return(pA)
}
doc_df = fread('../Data/Dissimilarity_Overlap.csv')
doc_df = doc_df[Dissimilarity!=0 & Same_Environment == TRUE & Carbon_Source_1 == 'Glucose' & Same_Inoculum==TRUE,]
p1<- plot_comp(doc_df,0.9,fn= '../Stat_Outputs/TTest_CSource_Cit_Overlap.csv')
p2<- plot_comp(doc_df,0.91,fn= '../Stat_Outputs/TTest_CSource_Cit_Overlap.csv')
p3<- plot_comp(doc_df,0.92,fn= '../Stat_Outputs/TTest_CSource_Cit_Overlap.csv')
p4<- plot_comp(doc_df,0.93,fn= '../Stat_Outputs/TTest_CSource_Cit_Overlap.csv')
p5<- plot_comp(doc_df,0.94,fn= '../Stat_Outputs/TTest_CSource_Cit_Overlap.csv')
p6<- plot_comp(doc_df,0.95,fn= '../Stat_Outputs/TTest_CSource_Cit_Overlap.csv')
p7<- plot_comp(doc_df,0.96,fn= '../Stat_Outputs/TTest_CSource_Cit_Overlap.csv')
p8<- plot_comp(doc_df,0.97,fn= '../Stat_Outputs/TTest_CSource_Cit_Overlap.csv')
p9<- plot_comp(doc_df,0.98,fn= '../Stat_Outputs/TTest_CSource_Cit_Overlap.csv')
p10<- plot_comp(doc_df,0.99,fn= '../Stat_Outputs/TTest_CSource_Cit_Overlap.csv')
p11<- plot_comp(doc_df,0.9925,fn= '../Stat_Outputs/TTest_CSource_Cit_Overlap.csv')
p12<- plot_comp(doc_df,0.995,fn= '../Stat_Outputs/TTest_CSource_Cit_Overlap.csv')

# p = grid.arrange(pA,pB,layout_matrix=rbind(c(1,1,1,1,1,2,2)))

ggsave('../Final_Figures/Supp6.png',grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol=3),height=12,width=9)
