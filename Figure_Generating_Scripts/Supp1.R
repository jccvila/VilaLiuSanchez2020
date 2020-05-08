library(data.table)
library(ggplot2)
library(ggpubr)
library(gridExtra)
rm(list=ls())

plot_lowess<- function(dat,fn,med){
  #Plot DOC using lowess curve fitting approach
  #med is the median overlap used,
  #fn is the file contained the boostrapped curve dta
  sdata= fread(fn)
  m = loess(dat$Dissimilarity ~ dat$Overlap,span=0.2,family="symmetric",iterations=5,statistics='none')
  sdata$y = predict(m,sdata$xs)
  sdata = sdata[!is.na(sdata$y)]
  p<- ggplot(dat) +
    geom_point(data=dat,aes(x= Overlap,y = Dissimilarity),size=1,alpha=0.05) +
    geom_segment(x=med,xend=med,y=0,yend=sqrt(log(2)),linetype=2,col='Red') +
    geom_ribbon(data=sdata,aes(x=xs,ymin=LCI_lowess,ymax=UCI_lowess),alpha=0.5,fill='Red') + 
    geom_line(data=sdata,aes(x=xs,y=y),col='Red',size=2) +
    scale_x_continuous(breaks=c(0,1),limits=c(0,1)) +
    scale_y_continuous(breaks=c(0,0.8),limits=c(0,0.83))  + theme(axis.line = element_line(size=1)) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.title=element_text(size=16))
  return(p)
}


plot_olm<- function(dat,fn,med){
  #Plot olm over subset of data with above med overlap
  # Fn is file containing bootstrapped data
  m = lm(Dissimilarity ~ Overlap, dat[Overlap>med])
  slope =floor(signif(m$coefficients[2],2)*100)/100
  sdata= fread(fn)
  dat = dat[Overlap>med,]
  sdata = sdata[xs>min(dat$Overlap),]
  pvalue = mean(sdata$PValue)
  if(pvalue ==0){pvalue = 'p < 0.002'}else{
    pvalue = paste('p = ',signif(pvalue,2),sep='')
  }
  
  p<- ggplot(dat) +
    geom_point(data=dat,aes(x= Overlap,y = Dissimilarity),size=1,alpha=0.05) +
    theme_pubr() + 
    geom_line(data=sdata,aes(x=xs,y=mean_lm),col='Blue',size=2) +
    geom_ribbon(data=sdata,aes(x=xs,ymin=LCI_lm,ymax=UCI_lm),alpha=0.5,fill='Blue') +
    scale_x_continuous(breaks=c(floor(med*100)/100,1),limits=c(floor(med*100)/100,1)) +
    scale_y_continuous(breaks=c(0,0.8),limits=c(0,0.9)) + theme(axis.line = element_line(size=1)) +
    annotate(geom='text',x = floor(min(dat$Overlap)*100)/100,y=0.9,label=paste('\t\t\t\t m =', slope,' ' , pvalue),size=2) +
    theme(plot.margin = margin(0, 0, 0, 0, "cm")) +
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  return(p)
}


plot_olm_cs<- function(dat,med){
  #Plot olm over subset of data with above med overlap gationg by carbon source for figure C and supplementary
  # Fn is file containing bootstrapped data
  #Subset
  dat = dat[Overlap>med]

  G = dat[Carbon_Source_1 == 'Glucose']
  C = dat[Carbon_Source_1 == 'Citrate']
  L = dat[Carbon_Source_1 == 'Leucine']
  #media
  mG = lm(Dissimilarity ~ Overlap, G)
  mC =lm(Dissimilarity ~ Overlap, C)
  mL =lm(Dissimilarity ~ Overlap, L)
  
  slopeG =floor(signif(mG$coefficients[2],2)*100)/100
  slopeC =floor(signif(mC$coefficients[2],2)*100)/100
  slopeL =floor(signif(mL$coefficients[2],2)*100)/100
  
  sdataG= fread( '../Stat_Outputs/Curve_Fitting_Glucose_Null.csv')
  sdataC= fread( '../Stat_Outputs/Curve_Fitting_Citrate_Null.csv')
  sdataL= fread( '../Stat_Outputs/Curve_Fitting_Leucine_Null.csv')
  
  sdataG = sdataG[xs>med,]
  pvalueG = mean(sdataG$PValue)
  sdataG$Carbon_Source_1 = 'Glucose'
  sdataC = sdataC[xs>med,]
  pvalueC = mean(sdataC$PValue)
  sdataC$Carbon_Source_1 = 'Citrate'
  sdataL = sdataL[xs>med,]
  sdataL$Carbon_Source_1 = 'Leucine'
  pvalueL = mean(sdataL$PValue)
  if(pvalueG ==0){pvalueG = 'p < 0.002'}else{
    pvalueG = paste('p = ',signif(pvalueG,2),sep='')
  }
  if(pvalueC ==0){pvalueC = 'p < 0.002'}else{
    pvalueC = paste('p = ',signif(pvalueC,2),sep='')
  }
  if(pvalueL ==0){pvalueL = 'p < 0.002'}else{
    pvalueL = paste('p = ',signif(pvalueL,2),sep='')
  }
  sdata = rbind(sdataG,sdataC,sdataL)
  lab_data = data.frame(Carbon_Source_1 = c('Glucose','Citrate','Leucine'),
                        lab = paste('\t\t\t\t m =', c(slopeG,slopeC,slopeL),' ' , c(pvalueG,pvalueC,pvalueL)) ,
                        Overlap = c(med,med,med),Dissimilarity = 0.9)
  dat$Carbon_Source_1 = factor(dat$Carbon_Source_1,levels=c('Glucose','Citrate','Leucine'))
  lab_data$Carbon_Source_1 = factor(lab_data$Carbon_Source_1,levels=c('Glucose','Citrate','Leucine'))
  p<- ggplot(dat) +
    geom_point(data=dat,aes(x= Overlap,y = Dissimilarity),size=1,alpha=0.05) +
    theme_bw() + 
    geom_line(data=sdata,aes(x=xs,y=mean_lm),col='Blue',size=2) +
    geom_ribbon(data=sdata,aes(x=xs,ymin=LCI_lm,ymax=UCI_lm),alpha=0.5,fill='Blue') +
    scale_y_continuous(breaks=c(0,0.8),limits=c(0,0.9)) +# theme(axis.line = element_line(size=1)) +
    geom_text(data = lab_data,size=2,aes(x=Overlap,y=Dissimilarity,label=lab)) +
    facet_wrap(~Carbon_Source_1,ncol=1)  +
    scale_x_continuous(limits=c(floor(min(dat$Overlap)*100)/100,1),breaks=c(floor(min(dat$Overlap)*100)/100,1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(p)
}

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
    stat_pvalue_manual(pvals,label='adj',y.position = c(0.85,0.9,0.95),size=3) +
    theme(legend.position = '' ,
          axis.line = element_line(size=1),
          axis.text = element_text(size=8)) +    
    scale_x_discrete(labels=c(expression('Glc-Glc'),expression('Cit-Cit'),expression('Glc-Cit')))
  return(p1)
}


doc_df = fread('../Data/Dissimilarity_Overlap_Null.csv')
doc_df = doc_df[Dissimilarity>0,]
doc_df2 = doc_df[Same_Environment ==TRUE,]

pB <- plot_lowess(doc_df2,'../Stat_Outputs/Curve_Fitting_SameEnv_Null.csv',median(doc_df2$Overlap))
pBInset <- plot_olm(doc_df2,'../Stat_Outputs/Curve_Fitting_SameEnv_Null.csv',median(doc_df2$Overlap))+  
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),axis.text =  element_blank(),axis.title = element_blank())
pC <- plot_olm_cs(doc_df2,median(doc_df2$Overlap))
pD <- plot_olm(doc_df2[Same_Inoculum==TRUE],'../Stat_Outputs/Curve_Fitting_SameEnvSameInoc_Null.csv',median(doc_df2$Overlap))
pE <- plot_olm(doc_df2[Same_Inoculum==FALSE],'../Stat_Outputs/Curve_Fitting_SameEnvDiffInoc_Null.csv',median(doc_df2$Overlap)) 
  
pF <- plot_comp_self(doc_df2,0.98)
pG <- plot_overlap_dists(doc_df)
pH <- plot_comp_diff(doc_df,0.98,fn = '../Stat_Outputs/TTest_CSource_Null.csv')

t1 = text_grob(expression('' == 'Inoculum'),face='bold',size=12)
t2 = text_grob(expression('' != 'Inoculum'),face='bold',size=12)
pleft = ggarrange(pB + annotation_custom(arrangeGrob(pBInset),xmin=0.0,xmax=0.3,ymin=0.6,ymax=0.83),labels='A')
pmiddle = ggarrange(pC,labels='B')
pright = ggarrange(grid.arrange(pD,top=t1),grid.arrange(pE,top=t2),ncol=1,labels=c('C','D'))
pbottom = ggarrange(pF,pG,pH,labels=c('E','F','G'),ncol=3)
ggsave(filename = '../Final_Figures/Supp1.png',grid.arrange(pleft,pmiddle,pright,pbottom,
                                                           layout_matrix = rbind(c(1,1,2,3),
                                                                                 c(1,1,2,3),
                                                                                 c(1,1,2,3),
                                                                                 c(4,4,4,4),
                                                                                 c(4,4,4,4))),height=8,width=10)