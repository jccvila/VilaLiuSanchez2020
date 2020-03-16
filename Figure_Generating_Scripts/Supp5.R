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

plot_comp <- function(doc_df,t){
  Ad <- doc_df[Carbon_Source_1 == 'Glucose'& Same_Environment == TRUE & Overlap>t,]
  Bd <- doc_df[Carbon_Source_1 == 'Citrate'& Same_Environment == TRUE & Overlap>t,]
  Cd <- doc_df[Carbon_Source_2 == 'Citrate'& Carbon_Source_1 =='Glucose' & Overlap>t,]
  A = rbind(Ad,Bd,Cd)
  A$Treatment = 'NA'
  A[Same_Environment == TRUE & Carbon_Source_1 == 'Glucose',]$Treatment = 'A'
  A[Same_Environment == TRUE & Carbon_Source_1 == 'Citrate',]$Treatment = 'B'
  A[Same_Environment == FALSE,]$Treatment ='C'
  mycomparisons = list( c('A','B'),c('B','C'),c('A','C'))
  
  p1 <- ggboxplot(A,x='Treatment',y='Dissimilarity',col='Treatment',palette = "jco",
                  add = "jitter",legend='right') +
    stat_compare_means(comparisons = mycomparisons,method='t.test') + scale_y_continuous(breaks=c(0,0.8)) +
    stat_compare_means(label.y = 1.2,method = "anova") + guides(col=FALSE) + 
    theme(axis.title.x=element_blank(),
          legend.position = '' ,
          axis.line = element_line(size=1),
          axis.text = element_text(size=8)) +    ggtitle(paste('t =', t)) +
    scale_x_discrete(labels=c(expression('Glc-Glc'),expression('Cit-Cit'),expression('Glc-Cit')))
  return(p1)
}

doc_df = fread('../Data/Dissimilarity_Overlap.csv')
doc_df = doc_df[Dissimilarity!=0,]
p1<- plot_comp(doc_df,0.9)
p2<- plot_comp(doc_df,0.91)
p3<- plot_comp(doc_df,0.92)
p4<- plot_comp(doc_df,0.93)
p5<- plot_comp(doc_df,0.94)
p6<- plot_comp(doc_df,0.95)
p7<- plot_comp(doc_df,0.96)
p8<- plot_comp(doc_df,0.97)
p9<- plot_comp(doc_df,0.98)
p10<- plot_comp(doc_df,0.99)
p11<- plot_comp(doc_df,0.995)
p12<- plot_comp(doc_df,0.999)

# p = grid.arrange(pA,pB,layout_matrix=rbind(c(1,1,1,1,1,2,2)))

ggsave('../Final Figures/Supp5.png',grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol=3),height=12,width=9)
