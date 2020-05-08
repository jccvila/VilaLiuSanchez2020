library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)
rm(list=ls())

plot_comp_self <- function(doc_df,t){
  Ad <- doc_df[Carbon_Source_1 == 'Glucose'& Overlap>t,]
  Bd <- doc_df[Carbon_Source_1 == 'Citrate' & Overlap>t,]
  Cd <- doc_df[Carbon_Source_1 == 'Leucine'& Overlap>t,]
  A = rbind(Ad,Bd,Cd)
  A$Treatment = A$Carbon_Source_1
  p1 <- ggboxplot(A,x='Treatment',y='Dissimilarity',col='Treatment',palette = "jco",
                  add = "jitter",legend='right') +
    labs(x = '') +
    theme(legend.position = '' ,
          axis.line = element_line(size=1),
          axis.text = element_text(size=8)) +    scale_y_continuous(breaks=c(0,0.8),limits=c(0,0.9))+
    scale_color_manual(values=c('#0073C2FF','#EFC000FF','Light Green')) +
    scale_x_discrete(labels=c(expression('Glucose'),expression('Citrate'),expression('Leucine'))) +
    geom_hline(yintercept=sqrt(log(2))/2,col='red',linetype=2)
  return(p1)
}


plot_overlap_dists <- function(doc_df){
  Sd = doc_df[Same_Environment==TRUE & Carbon_Source_1 %in% c('Glucose','Citrate'),]
  Dd = doc_df[Same_Environment==FALSE & Carbon_Source_1=='Glucose' & Carbon_Source_2 == 'Citrate',]
  p1 <- ggplot() + 
    geom_histogram(data=Dd,binwidth = .04,aes(x=Overlap,y=..count../sum(..count..)),fill='#868686FF') +
    geom_freqpoly(data=Sd,binwidth = .04,aes(x=Overlap,col=Carbon_Source_1,y=2*..count../sum(..count..)),linetype='dashed') +
    labs(x='O',y='P(O)') +guides(fill=FALSE,col=FALSE)+ 
    scale_color_manual(values=c('#EFC000FF','#0073C2FF')) +
    scale_x_continuous(breaks=c(0.0,1.0))  +
    scale_y_continuous(breaks=c(0.0,0.08)) +
    theme_classic() +
    theme(axis.text = element_text(size=8)) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(size=1))
  
  return(p1)
}

plot_comp_diff <- function(doc_df,o,fn){
  Ad <- doc_df[Carbon_Source_1 == 'Glucose' & Same_Environment == TRUE & Overlap>o,]
  Bd <- doc_df[Carbon_Source_1 == 'Citrate' &Same_Environment == TRUE & Overlap>o,]
  Cd <- doc_df[Carbon_Source_2 == 'Citrate'& Carbon_Source_1 =='Glucose' & Overlap>o,]
  A = rbind(Ad,Bd,Cd)
  A$Treatment = 'NA'
  A[Same_Environment == TRUE & Carbon_Source_1 == 'Glucose',]$Treatment = 'A'
  A[Same_Environment == TRUE & Carbon_Source_1 == 'Citrate',]$Treatment = 'B'
  A[Same_Environment == FALSE,]$Treatment ='C'
  
  tests = fread(fn)
  tests = tests[!is.na(tests$t) & Threshold==o,]
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
  pvals$adj = paste('p < ', 1/n)
  pvals[pvals$p !=0,]$adj = as.character(pvals[pvals$p !=0,]$p)
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
    scale_x_discrete(labels=c(expression('Glc-Glc'),expression('Cit-Cit'),expression('Glc-Cit')))
  return(p1)
}


doc_df = fread('../Data/Dissimilarity_Overlap.csv')
doc_df = doc_df[Dissimilarity>0& Same_Inoculum==TRUE,]
doc_df2 = doc_df[Same_Environment ==TRUE,]

pF <- plot_comp_self(doc_df2,0.98)
pG <- plot_overlap_dists(doc_df)
pH <- plot_comp_diff(doc_df,0.98,fn = '../Stat_Outputs/TTest_CSource_Same_Inoc.csv')
ggsave('../Final_Figures/Supp4.png',ggarrange(pF,pG,pH,labels='AUTO',nrow=1),height=3,width=10)

