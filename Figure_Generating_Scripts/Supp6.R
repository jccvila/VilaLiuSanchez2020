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
ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 
'%!in%' <- function(x,y)!('%in%'(x,y))

null = 0
if(null == 0){
  doc_df = fread('../Data/Dissimilarity_Overlap.csv')
}else{
  doc_df = fread('../Data/Dissimilarity_Overlap_Null.csv')
}

doc_df = doc_df[Carbon_Source_1 == 'Glucose' & Carbon_Source_2 == 'Glucose' & Same_Inoculum==TRUE]
# doc_df_b = doc_df_b[Carbon_Source_1 == 'Glucose' & Carbon_Source_2 == 'Glucose']
doc_df = doc_df[Dissimilarity!=0]
# doc_df_b = doc_df_b[Dissimilarity!=0]

plot_cit <- function(dat,med){
  dat = doc_df[Same_Inoculum==TRUE,]
  odata = fread('../Data/Emergent_Comunity_Data_Equilibrium.csv')[Carbon_Source=='Glucose' & Transfer==12]
  cit_coms = unique(odata[ESV_ID == 'Citrobacter',]$Sample_ID)
  dat$Citrobacter <- 'One'
  dat[Var1 %!in% cit_coms & Var2 %!in% cit_coms,]$Citrobacter= 'Neither'
  dat[Var1 %in% cit_coms & Var2 %in% cit_coms,]$Citrobacter= 'Both'
  dat$Citrobacter = factor(dat$Citrobacter,levels=c('Both','One','Neither'))
  dat = dat[Overlap>med,]

  mycomparisons = list( c('Both','One'),c('One','Neither'),c('Both','Neither'))

  pA <- ggboxplot(dat,x='Citrobacter',y='Dissimilarity',col='Citrobacter',palette = "jco",
                  add = "jitter",legend='right') + labs(x='') +
    stat_compare_means(comparisons = mycomparisons,method='t.test') + scale_x_discrete(labels=c('','','')) +
    scale_y_continuous(breaks=c(0,0.8)) +
    stat_compare_means(label.y = 1.2,method = "anova") +  
    theme(axis.text.x = element_blank(),legend.title = element_blank())+ theme_pubr() + 
    theme(axis.title= element_text(size=16),axis.line = element_line(size=2),axis.text = element_text(size=14)) +
    ggtitle(paste('t =', med))
  return(pA)
}


p1<- plot_cit(doc_df,0.9)
p2<- plot_cit(doc_df,0.91)
p3<- plot_cit(doc_df,0.92)
p4<- plot_cit(doc_df,0.93)
p5<- plot_cit(doc_df,0.94)
p6<- plot_cit(doc_df,0.95)
p7<- plot_cit(doc_df,0.96)
p8<- plot_cit(doc_df,0.97)
p9<- plot_cit(doc_df,0.98)
p10<- plot_cit(doc_df,0.99)
p11<- plot_cit(doc_df,0.9925)
p12<- plot_cit(doc_df,0.995)


# p = grid.arrange(pA,pB,layout_matrix=rbind(c(1,1,1,1,1,2,2)))

ggsave('../Final Figures/Supp6.png',ggarrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol=3,nrow=4,common.legend=TRUE),height=12,width=9)

