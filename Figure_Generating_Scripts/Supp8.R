rm(list=ls())
library(data.table)
library(ggplot2)
library(ggpubr)

plot_high_overlap <- function(dat,med){
  odata = fread('../Data/Emergent_Comunity_Data_Equilibrium.csv')[Carbon_Source=='Citrate' & Transfer==12]
  cit_coms = unique(odata[ESV_ID == 'Raoultella',]$Sample_ID)
  dat$Raoultella <- 'One'
  dat[Var1 %!in% cit_coms & Var2 %!in% cit_coms,]$Raoultella= 'Neither'
  dat[Var1 %in% cit_coms & Var2 %in% cit_coms,]$Raoultella= 'Both'
  dat$Raoultella = factor(dat$Raoultella,levels=c('Both','One','Neither'))
  med =  floor(median(dat$Overlap)*100)/100
  dat = dat[Overlap>med,]
  pA <- ggplot() + geom_point(dat,mapping = aes(x=Overlap,y=Dissimilarity,col=Raoultella)) +
    theme_pubr() + guides(fill=FALSE) + 
    # scale_x_continuous(breaks=c(floor(min(dat$Overlap)*1000)/1000,1),limits=c(floor(min(dat$Overlap)*1000)/1000,1)) + 
    scale_x_continuous(breaks=c(floor(med*100)/100,1),limits=c(floor(med*100)/100,1)) + 
    guides(alpha=FALSE,size=FALSE) +
    scale_y_continuous(breaks=c(0,0.8),limits=c(0,sqrt(log(2)))) +# stat_ellipse(type="norm",level=0.5) +
    theme(axis.title= element_text(size=16),axis.line = element_line(size=2),
          axis.text = element_text(size=14),legend.text = element_text(size=12)) + labs(col='Raoultella In') +
    scale_color_manual(values=c('#0073C2FF','#8EE5EE','#EFC000FF')) 
  return(pA)
}

plot_comp <- function(dat,med,fn){
  odata = fread('../Data/Emergent_Comunity_Data_Equilibrium.csv')[Carbon_Source=='Citrate' & Transfer==12]
  cit_coms = unique(odata[ESV_ID == 'Raoultella',]$Sample_ID)
  dat$Raoultella <- 'One'
  dat[Var1 %!in% cit_coms & Var2 %!in% cit_coms,]$Raoultella= 'Neither'
  dat[Var1 %in% cit_coms & Var2 %in% cit_coms,]$Raoultella= 'Both'
  dat$Raoultella = factor(dat$Raoultella,levels=c('Both','One','Neither'))
  dat = dat[Overlap>med,]
  tests = fread(fn)
  tests = tests[Threshold ==med & !is.na(tests$t)]
  pvals = data.frame()
  for(x in unique(tests$Comparison)){
    pvals = rbind(pvals,data.frame(Comparison = x,p = sum(tests[Comparison == x]$t<0)/nrow(tests[Comparison == x])))
  }
  pvals$adj = 'p < 0.002'
  pvals[pvals$p !=0,]$adj = as.character(floor(pvals[pvals$p !=0,]$p*1000)/1000)
  pvals[pvals$p!=0]
  pvals$group1 = c('Both','Both','One')
  pvals$group2 = c('One','Neither','Neither')
  pB <- ggboxplot(dat,x='Raoultella',y='Dissimilarity',col='Raoultella',palette = "jco",
                  add = "jitter",legend='right') +
    # stat_compare_means(comparisons = mycomparisons,method='t.test',size=2) + scale_y_continuous(breaks=c(0,0.8)) +
    guides(col=FALSE) +   labs(x = '') +
    stat_pvalue_manual(pvals,label='adj',y.position = c(0.8,0.9,0.85),size=3) +
    theme(legend.position = '' ,
          axis.line = element_line(size=1),
          axis.text = element_text(size=8))
  return(pB)
}



doc_df = fread('../Data/Dissimilarity_Overlap.csv')
doc_df = doc_df[Dissimilarity!=0 & Same_Environment == TRUE & Carbon_Source_1 == 'Citrate',]
  pA <- plot_high_overlap(doc_df[Same_Inoculum==TRUE],
                        med = median(doc_df[Same_Inoculum==TRUE]$Overlap))
pB <- plot_comp(doc_df[Same_Inoculum==TRUE],
                med = median(doc_df[Same_Inoculum==TRUE]$Overlap),
                fn = '../Stat_Outputs/TTest_CSource_Raoul.csv')
pAB = ggarrange(pA,pB,labels = 'AUTO')

ggsave('../Final_Figures/Supp8.png',pAB,width=9,height=4)