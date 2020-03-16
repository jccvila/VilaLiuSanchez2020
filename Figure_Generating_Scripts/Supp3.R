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

plot_heatmap<- function(data,On,Dn,med,spv,fpv){
  m = lm(Dissimilarity ~ Overlap, data[Overlap>med])
  slope =floor(signif(m$coefficients[2],2)*100)/100
  xs = seq(med,1,0.01)
  ys = xs*m$coefficients[2] + m$coefficients[1]
    sdata = data.frame(xs,ys)
  lr = sum(data$Overlap>0.5 & data$Dissimilarity<sqrt(log(2))/2)
  lr_null = c()
  for(j in 1:10000){
    lr_null = c(lr_null,sum(sample(On,nrow(data),replace=FALSE)>0.5 & sample(Dn,nrow(data),replace=FALSE)<sqrt(log(2))/2))
  }
  pvalue = length(which(lr_null>=lr))/length(lr_null)
  frac = signif(lr/nrow(data),2)
  if(fpv < 0.01 & frac>0.5){pvalue <- paste('f = ',frac,'* ',sep='')} else(pvalue <- paste('f = ',frac,' ,',sep=''))
  if(spv < 0.01){pvalue <- paste(pvalue,' m = ',slope,'*',sep='')} else(pvalue <- paste(pvalue,' slope = ',slope))
  p <- ggplot(data) + geom_point(data=data,aes(x=Overlap,y=Dissimilarity),alpha=0) +
    geom_bin2d(data = data, aes(x=Overlap, y=Dissimilarity),binwidth=c(0.05,0.05)) + 
    theme_pubr() + guides(fill=FALSE) + 
    scale_fill_gradient2(low='White',mid='Yellow',high='Red') +
    geom_line(data = sdata,aes(xs,ys),linetype=2) +
    scale_x_continuous(breaks = c(0,1),limits=c(-5.0e-09,10.5e-01),labels=c('','')) + 
    scale_y_continuous(breaks=c(0,0.8),limits=c(-5.0e-09,8.5e-01)) + 
    geom_text (x = 0.5,y = 0.86,label = pvalue,size=3)
  return(p)
}

plot_overlap_dists <- function(doc_df){
  Sd = doc_df[Same_Environment==TRUE,]
  Sd$Carbon_Source_1 = factor(Sd$Carbon_Source_1,levels=c('Glucose','Citrate','Leucine'))
  p1 <- ggplot(data=Sd, aes(x=Overlap,fill=Carbon_Source_1)) + 
    geom_histogram(binwidth = .04,aes(y=..count../sum(..count..))) +
    labs(x='O',y='P(O)') +guides(fill=FALSE)+
    scale_fill_manual(values=c('#0073C2FF','#EFC000FF','Light Green')) +
    theme_classic() + 
    facet_wrap(~Carbon_Source_1,nrow=3,scales='fixed') + theme_pubr() + scale_x_continuous(breaks=c(0.0,1.0)) +
    scale_y_continuous(breaks=c(0.0,0.05)) +    
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
    geom_histogram(data=Dd,binwidth = .04,aes(x=Overlap,fill=Treatment,y=..count../sum(..count..))) +
    geom_freqpoly(data=Sd,binwidth = .04,aes(x=Overlap,col=Carbon_Source_1,y=2*..count../sum(..count..)),linetype='dashed') +
    labs(x='O',y='') +guides(fill=FALSE,col=FALSE)+ 
    scale_fill_manual(values=c('#868686FF','#A73030FF','#20854EFF')) +
    scale_color_manual(values=c('#0073C2FF','#EFC000FF','Light Green')) +
    geom_hline(yintercept=0, colour="white", size=1) +
    facet_wrap(~Treatment,nrow=3,scales='fixed') + theme_pubr() + scale_x_continuous(breaks=c(0.0,1.0))  +
    scale_y_continuous(breaks=c(0.0,0.09)) +    
    theme(axis.text = element_text(size=8))

  return(ggarrange(p1,p2,labels='AUTO'))
}

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
          axis.text = element_text(size=8)) +    
    scale_x_discrete(labels=c(expression('Glc-Glc'),expression('Cit-Cit'),expression('Glc-Cit')))
  return(ggarrange(p1,labels=c('c')))
}

add_vp <- function(x){return(grid.arrange(x,vp=viewport(width=0.95, height=0.95)))}
null = 0
doc_df = fread('../Data/Dissimilarity_Overlap.csv')
doc_df[Dissimilarity>0,]
pAB <- plot_overlap_dists(doc_df)
pC<- plot_comp(doc_df,0.98)
# p = grid.arrange(pA,pB,layout_matrix=rbind(c(1,1,1,1,1,2,2)))

ggsave('../Final Figures/Supp3.png',grid.arrange(pAB,pC,ncol=2),height=4.5,width=9)
