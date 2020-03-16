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
library(cowplot)
set.seed(1)

g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

plot_overlap_dists <- function(doc_df){
  Sd = doc_df[Same_Environment==TRUE,]
  Sd$Carbon_Source_1 = factor(Sd$Carbon_Source_1,levels=c('Glucose','Citrate','Leucine'))
  p1 <- ggplot(data=Sd, aes(x=Overlap,fill=Carbon_Source_1)) + 
    geom_histogram(binwidth = .04,aes(y=..count../sum(..count..))) +
    labs(x='O',y='P(O)') +guides(fill=FALSE)+
    scale_fill_manual(values=c('#0073C2FF','#EFC000FF','Light Green')) +
    geom_bar(aes(y=..prop.., group = 1)) + theme_classic() + 
    facet_wrap(~Carbon_Source_1,nrow=3,scales='fixed') + theme_pubr() + scale_x_continuous(breaks=c(0.0,1.0)) +
    scale_y_continuous(breaks=c(0.0,0.15)) +    
    theme(axis.text = element_text(size=8))
  
  
  Dd = doc_df[Same_Environment==FALSE,]
  Dd[,Treatment:= '']
  Dd[Carbon_Source_1=='Glucose' & Carbon_Source_2 == 'Citrate',]$Treatment = 'Glc-Cit'
  Dd[Carbon_Source_1=='Leucine' & Carbon_Source_2 == 'Citrate',]$Treatment = 'Cit-Leu'
  Dd[Carbon_Source_1=='Leucine' & Carbon_Source_2 == 'Glucose',]$Treatment = 'Glc-Leu'
  Dd$Treatment = factor(Dd$Treatment,levels=c('Glc-Cit','Glc-Leu','Cit-Leu'))
  Sd[,Treatment:= '']
  Sd2 = Sd
  Sd[Carbon_Source_1=='Glucose',]$Treatment = 'Glc-Cit'
  Sd[Carbon_Source_1=='Citrate',]$Treatment = 'Cit-Leu'
  Sd[Carbon_Source_1=='Leucine',]$Treatment = 'Glc-Leu' 
  Sd2[Carbon_Source_1=='Citrate',]$Treatment = 'Glc-Cit'
  Sd2[Carbon_Source_1=='Leucine',]$Treatment = 'Cit-Leu'
  Sd2[Carbon_Source_1=='Glucose',]$Treatment = 'Glc-Leu'
  Sd$Treatment = factor(Sd$Treatment,levels=c('Glc-Cit','Glc-Leu','Cit-Leu'))
  
  Sd = rbind(Sd,Sd2)
  p2 <- ggplot() + 
    geom_histogram(data=Dd[Treatment=='Glc-Cit'],binwidth = .04,aes(x=Overlap,fill=Treatment,y=..count../sum(..count..))) +
    geom_freqpoly(data=Sd[Treatment=='Glc-Cit'],binwidth = .04,aes(x=Overlap,col=Carbon_Source_1,y=2*..count../sum(..count..)),linetype='dashed') +
    labs(x='O',y='P(O)') +guides(fill=FALSE,col=FALSE)+ 
    scale_fill_manual(values=c('#868686FF','#A73030FF','#20854EFF')) +
    scale_color_manual(values=c('#0073C2FF','#EFC000FF','Light Green')) +
    geom_hline(yintercept=0, colour="white", size=1) +
    #facet_wrap(~Treatment,nrow=3,scales='fixed') + #theme_bw() +
    scale_x_continuous(breaks=c(0.0,1.0))  +
    scale_y_continuous(breaks=c(0.0,0.5)) +theme_classic() +
    theme(axis.text = element_text(size=8)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(size=1))
  
  return(ggarrange(p2,labels='B'))
}


plot_overlap_dists2 <- function(doc_df){
  Sd = doc_df[Same_Environment==TRUE & Same_Inoculum==TRUE,]
  Sd$Carbon_Source_1 = factor(Sd$Carbon_Source_1,levels=c('Glucose','Citrate','Leucine'))
  p1 <- ggplot(data=Sd, aes(x=Overlap,fill=Carbon_Source_1)) + 
    geom_histogram(binwidth = .04,aes(y=..count../sum(..count..))) +
    labs(x='O',y='P(O)') +guides(fill=FALSE)+
    scale_fill_manual(values=c('#0073C2FF','#EFC000FF','Light Green')) +
    geom_bar(aes(y=..prop.., group = 1)) + theme_classic() + 
    facet_wrap(~Carbon_Source_1,nrow=3,scales='fixed') + theme_pubr() + scale_x_continuous(breaks=c(0.0,1.0)) +
    scale_y_continuous(breaks=c(0.0,0.15)) +    
    theme(axis.text = element_text(size=8))
  
  
  Dd = doc_df[Same_Environment==FALSE  & Same_Inoculum==TRUE,]
  Dd[,Treatment:= '']
  Dd[Carbon_Source_1=='Glucose' & Carbon_Source_2 == 'Citrate',]$Treatment = 'Glc-Cit'
  Dd[Carbon_Source_1=='Leucine' & Carbon_Source_2 == 'Citrate',]$Treatment = 'Cit-Leu'
  Dd[Carbon_Source_1=='Leucine' & Carbon_Source_2 == 'Glucose',]$Treatment = 'Glc-Leu'
  Dd$Treatment = factor(Dd$Treatment,levels=c('Glc-Cit','Glc-Leu','Cit-Leu'))
  Sd[,Treatment:= '']
  Sd2 = Sd
  Sd[Carbon_Source_1=='Glucose',]$Treatment = 'Glc-Cit'
  Sd[Carbon_Source_1=='Citrate',]$Treatment = 'Cit-Leu'
  Sd[Carbon_Source_1=='Leucine',]$Treatment = 'Glc-Leu' 
  Sd2[Carbon_Source_1=='Citrate',]$Treatment = 'Glc-Cit'
  Sd2[Carbon_Source_1=='Leucine',]$Treatment = 'Cit-Leu'
  Sd2[Carbon_Source_1=='Glucose',]$Treatment = 'Glc-Leu'
  Sd$Treatment = factor(Sd$Treatment,levels=c('Glc-Cit','Glc-Leu','Cit-Leu'))
  
  Sd = rbind(Sd,Sd2)
  p2 <- ggplot() + 
    geom_histogram(data=Dd[Treatment=='Glc-Cit'],binwidth = .04,aes(x=Overlap,fill=Treatment,y=..count../sum(..count..))) +
    geom_freqpoly(data=Sd[Treatment=='Glc-Cit'],binwidth = .04,aes(x=Overlap,col=Carbon_Source_1,y=2*..count../sum(..count..)),linetype='dashed') +
    labs(x='O',y='P(O)') +guides(fill=FALSE,col=FALSE)+ 
    scale_fill_manual(values=c('#868686FF','#A73030FF','#20854EFF')) +
    scale_color_manual(values=c('#0073C2FF','#EFC000FF','Light Green')) +
    geom_hline(yintercept=0, colour="white", size=1) +
    #facet_wrap(~Treatment,nrow=3,scales='fixed') + #theme_bw() +
    scale_x_continuous(breaks=c(0.0,1.0))  +
    scale_y_continuous(breaks=c(0.0,0.5)) +theme_classic() +
    theme(axis.text = element_text(size=8)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(size=1))
  
  return(ggarrange(p2,labels='B'))
}


plot_comp_self <- function(doc_df,t){
  Ad <- doc_df[Carbon_Source_1 == 'Glucose'& Same_Environment == TRUE & Overlap>t,]
  Bd <- doc_df[Carbon_Source_1 == 'Citrate'& Same_Environment == TRUE & Overlap>t,]
  Cd <- doc_df[Carbon_Source_1 == 'Leucine'& Same_Environment ==TRUE & Overlap>t,]
  A = rbind(Ad,Bd,Cd)
  A$Treatment = A$Carbon_Source_1
  p1 <- ggboxplot(A,x='Treatment',y='Dissimilarity',col='Treatment',palette = "jco",
                  add = "jitter",legend='right') +
    # stat_compare_means(comparisons = mycomparisons,method='t.test') + scale_y_continuous(breaks=c(0,0.8)) +
    # stat_compare_means(label.y = 1.2,method = "anova") + guides(col=FALSE) + 
    labs(x = '') +
    theme(legend.position = '' ,
          axis.line = element_line(size=1),
          axis.text = element_text(size=8)) +    scale_y_continuous(breaks=c(0,0.8),limits=c(0,0.9))+
    scale_color_manual(values=c('#0073C2FF','#EFC000FF','Light Green')) +
    scale_x_discrete(labels=c(expression('Glucose'),expression('Citrate'),expression('Leucine'))) +
    geom_hline(yintercept=sqrt(log(2))/2,col='red',linetype=2)
  # B = A[Same_Inoculum==TRUE,]
  # p2 <- ggboxplot(B,x='Treatment',y='Dissimilarity',col='Treatment',palette = "jco",
  #                 add = "jitter",legend='right') +
  #   # stat_compare_means(comparisons = mycomparisons,method='t.test') + scale_y_continuous(breaks=c(0,0.8)) +
  #   # stat_compare_means(label.y = 1.2,method = "anova") + guides(col=FALSE) + 
  #   scale_color_manual(values=c('#0073C2FF','#EFC000FF','Light Green')) +
  #   theme(axis.title.x=element_blank(),
  #         legend.position = '' ,
  #         axis.line = element_line(size=1),
  #         axis.text = element_text(size=8))  +    scale_y_continuous(breaks=c(0,0.8),limits=c(0,0.9))+
  #   scale_x_discrete(labels=c(expression('Glucose'),expression('Citrate'),expression('Leucine'))) +geom_hline(yintercept=sqrt(log(2))/2,col='red',linetype=2) 
  table(paste(A$Treatment,A$Dissimilarity>sqrt(log(2))/2))
  # table(paste(B$Treatment,B$Dissimilarity>sqrt(log(2))/2))
  return(ggarrange(p1,labels='A',ncol=1))
}


plot_comp_self2 <- function(doc_df,t){
  Ad <- doc_df[Carbon_Source_1 == 'Glucose'& Same_Environment == TRUE & Overlap>t,]
  Bd <- doc_df[Carbon_Source_1 == 'Citrate'& Same_Environment == TRUE & Overlap>t,]
  Cd <- doc_df[Carbon_Source_1 == 'Leucine'& Same_Environment ==TRUE & Overlap>t,]
  A = rbind(Ad,Bd,Cd)
  A$Treatment = A$Carbon_Source_1
  p1 <- ggboxplot(A,x='Treatment',y='Dissimilarity',col='Treatment',palette = "jco",
                  add = "jitter",legend='right') +
    # stat_compare_means(comparisons = mycomparisons,method='t.test') + scale_y_continuous(breaks=c(0,0.8)) +
    # stat_compare_means(label.y = 1.2,method = "anova") + guides(col=FALSE) + 
    labs(x = '') +
    theme(legend.position = '' ,
          axis.line = element_line(size=1),
          axis.text = element_text(size=8)) +    scale_y_continuous(breaks=c(0,0.8),limits=c(0,0.9))+
    scale_color_manual(values=c('#0073C2FF','#EFC000FF','Light Green')) +
    scale_x_discrete(labels=c(expression('Glucose'),expression('Citrate'),expression('Leucine'))) +
    geom_hline(yintercept=sqrt(log(2))/2,col='red',linetype=2)
  B = A[Same_Inoculum==TRUE,]
  p2 <- ggboxplot(B,x='Treatment',y='Dissimilarity',col='Treatment',palette = "jco",
                  add = "jitter",legend='right') +
    # stat_compare_means(comparisons = mycomparisons,method='t.test') + scale_y_continuous(breaks=c(0,0.8)) +
    # stat_compare_means(label.y = 1.2,method = "anova") + guides(col=FALSE) + 
    scale_color_manual(values=c('#0073C2FF','#EFC000FF','Light Green')) + labs(x='') +
    theme(legend.position = '' ,
          axis.line = element_line(size=1),
          axis.text = element_text(size=8))  +    scale_y_continuous(breaks=c(0,0.8),limits=c(0,0.9))+
    scale_x_discrete(labels=c(expression('Glucose'),expression('Citrate'),expression('Leucine'))) +geom_hline(yintercept=sqrt(log(2))/2,col='red',linetype=2) 
  table(paste(A$Treatment,A$Dissimilarity>sqrt(log(2))/2))
  table(paste(B$Treatment,B$Dissimilarity>sqrt(log(2))/2))
  return(ggarrange(p2,labels='A',ncol=1))
}


plot_comp_diff <- function(doc_df,t){
  Ad <- doc_df[Carbon_Source_1 == 'Glucose' & Same_Environment == TRUE & Overlap>t,]
  Bd <- doc_df[Carbon_Source_1 == 'Citrate' &Same_Environment == TRUE & Overlap>t,]
  Cd <- doc_df[Carbon_Source_2 == 'Citrate'& Carbon_Source_1 =='Glucose' & Overlap>t,]
  A = rbind(Ad,Bd,Cd)
  A$Treatment = 'NA'
  A[Same_Environment == TRUE & Carbon_Source_1 == 'Glucose',]$Treatment = 'A'
  A[Same_Environment == TRUE & Carbon_Source_1 == 'Citrate',]$Treatment = 'B'
  A[Same_Environment == FALSE,]$Treatment ='C'
  mycomparisons = list( c('A','B'),c('B','C'),c('A','C'))
  
  p1 <- ggboxplot(A,x='Treatment',y='Dissimilarity',col='Treatment',palette = "jco",
                  add = "jitter",legend='right') +
    stat_compare_means(comparisons = mycomparisons,method='t.test',size=2) + scale_y_continuous(breaks=c(0,0.8)) +
    stat_compare_means(label.y = 1.2,method = "anova",size=2) + guides(col=FALSE) +   labs(x = '') +

    theme(legend.position = '' ,
          axis.line = element_line(size=1),
          axis.text = element_text(size=8)) +    
    scale_x_discrete(labels=c(expression('Glc-Glc'),expression('Cit-Cit'),expression('Glc-Cit')))
  return(ggarrange(p1,labels=c('C')))
}

plot_comp_diff2 <- function(doc_df,t){
  Ad <- doc_df[Carbon_Source_1 == 'Glucose' & Same_Environment == TRUE & Overlap>t,]
  Bd <- doc_df[Carbon_Source_1 == 'Citrate' &Same_Environment == TRUE  & Overlap>t,]
  Cd <- doc_df[Carbon_Source_2 == 'Citrate'& Carbon_Source_1 =='Glucose' & Overlap>t,]
  A = rbind(Ad,Bd,Cd)
  A$Treatment = 'NA'
  A[Same_Environment == TRUE & Carbon_Source_1 == 'Glucose',]$Treatment = 'A'
  A[Same_Environment == TRUE & Carbon_Source_1 == 'Citrate',]$Treatment = 'B'
  A[Same_Environment == FALSE,]$Treatment ='C'
  mycomparisons = list( c('A','B'),c('B','C'),c('A','C'))
  
  p1 <- ggboxplot(A,x='Treatment',y='Dissimilarity',col='Treatment',palette = "jco",
                  add = "jitter",legend='right') +
    stat_compare_means(comparisons = mycomparisons,method='t.test',size=2) + scale_y_continuous(breaks=c(0,0.8)) +
    stat_compare_means(label.y = 1.2,method = "anova",size=2) + guides(col=FALSE) +   labs(x = '') +
    
    theme(legend.position = '' ,
          axis.line = element_line(size=1),
          axis.text = element_text(size=8)) +    
    scale_x_discrete(labels=c(expression('Glc-Glc'),expression('Cit-Cit'),expression('Glc-Cit')))
  return(ggarrange(p1,labels=c('C')))
}

null = 0

doc_df = fread('../Data/Dissimilarity_Overlap.csv')
doc_df = doc_df[Dissimilarity!=0,]

# fig_3 = grid.arrange(plot_comp_self(doc_df,0.98),plot_overlap_dists(doc_df),plot_comp_diff(doc_df,0.98),ncol=3)
s4 = grid.arrange(plot_comp_self(doc_df[Same_Inoculum==TRUE,],0.98),plot_overlap_dists(doc_df[Same_Inoculum==TRUE,]),plot_comp_diff(doc_df[Same_Inoculum==TRUE,],0.98),ncol=3)
# ggsave('../Final Figures/Fig2.png',fig_3,height=2,width=6)
ggsave('../Final Figures/Supp4.png',s4,height=2,width=6)

# 
