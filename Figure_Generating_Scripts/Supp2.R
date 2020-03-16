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

plot_lowess<- function(dat,fn,med){
  sdata= fread(fn)
  m = loess(dat$Dissimilarity ~ dat$Overlap,span=0.2,family="symmetric",iterations=5,statistics='none')
  sdata$y = predict(m,sdata$xs)
  sdata = sdata[!is.na(sdata$y)]
  p<- ggplot(dat) +
    geom_point(data=dat,aes(x= Overlap,y = Dissimilarity),size=1,alpha=0.05) +
    theme_pubr()  +
    geom_segment(x=med,xend=med,y=0,yend=sqrt(log(2)),linetype=2,col='Red') +
    geom_ribbon(data=sdata,aes(x=xs,ymin=LCI_lowess,ymax=UCI_lowess),alpha=0.5,fill='Red') + 
    geom_line(data=sdata,aes(x=xs,y=y),col='Red',size=2) +
    scale_x_continuous(breaks=c(0,1),limits=c(0,1)) +
    scale_y_continuous(breaks=c(0,0.8),limits=c(0,0.83))  + theme(axis.line = element_line(size=1))
  return(p)
}

plot_olm<- function(dat,fn,med){
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
    annotate(geom='text',x = floor(min(dat$Overlap)*100)/100,y=0.9,label=paste('\t\t\t\t m =', slope,' ' , pvalue),size=2)
  return(p)
}




plot_olm2<- function(dat,null){

  G = dat[Carbon_Source_1 == 'Glucose']
  C = dat[Carbon_Source_1 == 'Citrate']
  L = dat[Carbon_Source_1 == 'Leucine']
  
  med_G = median(dat$Overlap)
  med_C = median(dat$Overlap)
  med_L = median(dat$Overlap)
  
  mG = lm(Dissimilarity ~ Overlap, G[Overlap>median(med_G)])
  mC =lm(Dissimilarity ~ Overlap, C[Overlap>median(med_C)])
  mL =lm(Dissimilarity ~ Overlap, L[Overlap>median(med_L)])
  
  slopeG =floor(signif(mG$coefficients[2],2)*100)/100
  slopeC =floor(signif(mC$coefficients[2],2)*100)/100
  slopeL =floor(signif(mL$coefficients[2],2)*100)/100
  if(null==0){
    sdataG= fread( '../Stat_Outputs/DOC_stats_A.csv')
    sdataC= fread( '../Stat_Outputs/DOC_stats_B.csv')
    sdataL= fread( '../Stat_Outputs/DOC_stats_C.csv')
  }else{
    sdataG= fread( '../Stat_Outputs/DOC_stats_A_Null.csv')
    sdataC= fread( '../Stat_Outputs/DOC_stats_B_Null.csv')
    sdataL= fread( '../Stat_Outputs/DOC_stats_C_Null.csv')
  }
  sdataG = sdataG[xs>min(dat$Overlap),]
  pvalueG = mean(sdataG$PValue)
  sdataG$Carbon_Source_1 = 'Glucose'
  sdataC = sdataC[xs>min(dat$Overlap),]
  pvalueC = mean(sdataC$PValue)
  sdataC$Carbon_Source_1 = 'Citrate'
  sdataL = sdataL[xs>min(dat$Overlap),]
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
  dat$Carbon_Source_1 = factor(dat$Carbon_Source_1,levels=c('Glucose','Citrate','Leucine'))
  sdata = rbind(sdataG,sdataC,sdataL)
  lab_data = data.frame(Carbon_Source_1 = c('Glucose','Citrate','Leucine'),
                        lab = paste('\t\t\t\t m =', c(slopeG,slopeC,slopeL),' ' , c(pvalueG,pvalueC,pvalueL)) ,
                        Overlap = c(med_G,med_C,med_L),Dissimilarity = 0.9)
  dat$Carbon_Source_1 = factor(dat$Carbon_Source_1,levels=c('Glucose','Citrate','Leucine'))
  lab_data$Carbon_Source_1 = factor(lab_data$Carbon_Source_1,levels=c('Glucose','Citrate','Leucine'))
  sdata$Carbon_Source_1 = factor(sdata$Carbon_Source_1,levels=c('Glucose','Citrate','Leucine'))
  tdat = dat[(Carbon_Source_1 == 'Glucose' & Overlap>med_G) | (Carbon_Source_1 == 'Citrate' & Overlap>med_C)|  (Carbon_Source_1 == 'Leucine' & Overlap>med_L),]
  tsdata = sdata[(Carbon_Source_1 == 'Glucose' & xs>med_G) | (Carbon_Source_1 == 'Citrate' & xs>med_C)|  (Carbon_Source_1 == 'Leucine' & xs>med_L),]
  
  p<- ggplot(tdat) +
    geom_point(data=tdat,aes(x= Overlap,y = Dissimilarity),size=1,alpha=0.05) +
    theme_bw() + 
    geom_line(data=tsdata,aes(x=xs,y=mean_lm),col='Blue',size=2) +
    geom_ribbon(data=tsdata,aes(x=xs,ymin=LCI_lm,ymax=UCI_lm),alpha=0.5,fill='Blue') +
    scale_y_continuous(breaks=c(0,0.8),limits=c(0,0.9)) +# theme(axis.line = element_line(size=1)) +
    geom_text(data = lab_data,size=2,aes(x=Overlap,y=Dissimilarity,label=lab)) +
    facet_wrap(~Carbon_Source_1,ncol=1)  + scale_x_continuous(limits=c(0.54,1),breaks=c(0.54,1))
  return(p)
}

plot_olm2_2<- function(dat,null){
  
  G = dat[Carbon_Source_1 == 'Glucose' & Carbon_Source_2 == 'Citrate']
  C = dat[Carbon_Source_1 == 'Leucine' & Carbon_Source_2 == 'Citrate']
  L = dat[Carbon_Source_1 == 'Leucine' & Carbon_Source_2 == 'Glucose']
  dat$Carbon_Source = 'NA'
  dat[Carbon_Source_1 == 'Glucose' & Carbon_Source_2 == 'Citrate',]$Carbon_Source = 'Glc-Cit'
  dat[Carbon_Source_1 == 'Leucine' & Carbon_Source_2 == 'Citrate',]$Carbon_Source = 'Cit-Leu'
  dat[Carbon_Source_1 == 'Leucine' & Carbon_Source_2 == 'Glucose',]$Carbon_Source = 'Glc-Leu'
  med_G = median(dat$Overlap)
  med_C = median(dat$Overlap)
  med_L = median(dat$Overlap)
  
  mG = lm(Dissimilarity ~ Overlap, G[Overlap>median(med_G)])
  mC =lm(Dissimilarity ~ Overlap, C[Overlap>median(med_C)])
  mL =lm(Dissimilarity ~ Overlap, L[Overlap>median(med_L)])
  
  slopeG =floor(signif(mG$coefficients[2],2)*100)/100
  slopeC =floor(signif(mC$coefficients[2],2)*100)/100
  slopeL =floor(signif(mL$coefficients[2],2)*100)/100
  if(null==0){
    sdataG= fread( '../Stat_Outputs/DOC_stats_E.csv')
    sdataC= fread( '../Stat_Outputs/DOC_stats_F.csv')
    sdataL= fread( '../Stat_Outputs/DOC_stats_G.csv')
  }else{
    sdataG= fread( '../Stat_Outputs/DOC_stats_E_Null.csv')
    sdataC= fread( '../Stat_Outputs/DOC_stats_F_Null.csv')
    sdataL= fread( '../Stat_Outputs/DOC_stats_G_Null.csv')
  }
  sdataG = sdataG[xs>min(dat$Overlap),]
  pvalueG = mean(sdataG$PValue)
  sdataG$Carbon_Source = 'Glc-Cit'
  sdataC = sdataC[xs>min(dat$Overlap),]
  pvalueC = mean(sdataC$PValue)
  sdataC$Carbon_Source = 'Cit-Leu'
  sdataL = sdataL[xs>min(dat$Overlap),]
  sdataL$Carbon_Source = 'Glc-Leu'
  
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
  dat$Carbon_Source = factor(dat$Carbon_Source,levels=c('Glc-Cit','Cit-Leu','Glc-Leu'))
  sdata = rbind(sdataG,sdataC,sdataL)
  lab_data = data.frame(Carbon_Source = c('Glc-Cit','Cit-Leu','Glc-Leu'),
                        lab = paste('\t\t\t\t m =', c(slopeG,slopeC,slopeL),' ' , c(pvalueG,pvalueC,pvalueL)) ,
                        Overlap = c(med_G,med_C,med_L),Dissimilarity = 0.9)
  dat$Carbon_Source = factor(dat$Carbon_Source,levels=c('Glc-Cit','Cit-Leu','Glc-Leu'))
  lab_data$Carbon_Source = factor(lab_data$Carbon_Source,levels=c('Glc-Cit','Cit-Leu','Glc-Leu'))
  sdata$Carbon_Source = factor(sdata$Carbon_Source,levels=c('Glc-Cit','Cit-Leu','Glc-Leu'))
  tdat = dat[(Carbon_Source == 'Glc-Cit' & Overlap>med_G) | (Carbon_Source == 'Cit-Leu' & Overlap>med_C)|  (Carbon_Source == 'Glc-Leu' & Overlap>med_L),]
  tsdata = sdata[(Carbon_Source == 'Glc-Cit' & xs>med_G) | (Carbon_Source == 'Cit-Leu' & xs>med_C)|  (Carbon_Source == 'Glc-Leu' & xs>med_L),]
  
  p<- ggplot(tdat) +
    geom_point(data=tdat,aes(x= Overlap,y = Dissimilarity),size=1,alpha=0.05) +
    theme_bw() + 
    geom_line(data=tsdata,aes(x=xs,y=mean_lm),col='Blue',size=2) + 
    geom_ribbon(data=tsdata,aes(x=xs,ymin=LCI_lm,ymax=UCI_lm),alpha=0.5,fill='Blue') +
    scale_y_continuous(breaks=c(0,0.8),limits=c(0,0.9)) +# theme(axis.line = element_line(size=1)) +
    geom_text(data = lab_data,size=2,aes(x=Overlap,y=Dissimilarity,label=lab)) +
    facet_wrap(~Carbon_Source,ncol=1)  + scale_x_continuous(limits=c(0.42,1),breaks=c(0.42,1))
  return(p)
}

add_vp <- function(x){return(grid.arrange(x,vp=viewport(width=0.95, height=0.95)))}
null = 0

if(null == 0){
  doc_df = fread('../Data/Dissimilarity_Overlap.csv')
}else{
  doc_df = fread('../Data/Dissimilarity_Overlap_Null.csv')
}
doc_df = doc_df[Dissimilarity>0,]
b = text_grob(expression('Overlap'),face='bold',size=16)
r = text_grob(expression('Dissimilarity'),face='bold',size=16,rot=90)
doc_df2 = doc_df[Same_Environment ==TRUE,]
doc_df3 = doc_df[Same_Environment ==FALSE,]

med = median(doc_df2$Overlap)
med2 = median(doc_df3$Overlap)
if(null==0){
p1 <- plot_lowess(doc_df2,'../Stat_Outputs/DOC_stats.csv',med) +  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title=element_text(size=16))# +  #labs(x='',y='') +
  #theme(plot.margin = margin(0, 0, 0, 0, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1b <- plot_olm(doc_df2,'../Stat_Outputs/DOC_stats.csv',med)+  theme_bw()  + 
  theme(axis.title=element_blank(),axis.ticks=element_blank(),axis.text = element_blank()) + # + #labs(x='',y='') +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),axis.text =  element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())


p1_2 <- plot_lowess(doc_df3,'../Stat_Outputs/DOC_stats_D.csv',med2) +  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title=element_text(size=16))# +  #labs(x='',y='') +
#theme(plot.margin = margin(0, 0, 0, 0, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1b_2 <- plot_olm(doc_df3,'../Stat_Outputs/DOC_stats_D.csv',med2)+  theme_bw()  + 
  theme(axis.title=element_blank(),axis.ticks=element_blank(),axis.text = element_blank()) + # + #labs(x='',y='') +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),axis.text =  element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p3 <- plot_olm(doc_df2[Same_Inoculum==TRUE],'../Stat_Outputs/DOC_stats_H.csv',med)  + 
  theme(axis.title=element_text(size=16)) + theme_bw() + #labs(x='',y='') +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p4 <- plot_olm(doc_df2[Same_Inoculum==FALSE],'../Stat_Outputs/DOC_stats_I.csv',med)  + 
  theme(axis.title=element_text(size=16)) + theme_bw() + #labs(x='',y='') +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())


p3_2 <- plot_olm(doc_df3[Same_Inoculum==TRUE],'../Stat_Outputs/DOC_stats_J.csv',med2)  + 
  theme(axis.title=element_text(size=16)) + theme_bw() + #labs(x='',y='') +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

p4_2 <- plot_olm(doc_df3[Same_Inoculum==FALSE],'../Stat_Outputs/DOC_stats_K.csv',med2)  + 
  theme(axis.title=element_text(size=16)) + theme_bw() + #labs(x='',y='') +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}else{
  p1 <- plot_lowess(doc_df2,'../Stat_Outputs/DOC_stats_Null.csv',med) +  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.title=element_text(size=16))# +  #labs(x='',y='') +
  #theme(plot.margin = margin(0, 0, 0, 0, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p1b <- plot_olm(doc_df2,'../Stat_Outputs/DOC_stats_Null.csv',med)+  theme_bw()  + 
    theme(axis.title=element_blank(),axis.ticks=element_blank(),axis.text = element_blank()) + # + #labs(x='',y='') +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"),axis.text =  element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
  p1_2 <- plot_lowess(doc_df3,'../Stat_Outputs/DOC_stats_D_Null.csv',med2) +  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.title=element_text(size=16))# +  #labs(x='',y='') +
  #theme(plot.margin = margin(0, 0, 0, 0, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p1b_2 <- plot_olm(doc_df3,'../Stat_Outputs/DOC_stats_D_Null.csv',med2)+  theme_bw()  + 
    theme(axis.title=element_blank(),axis.ticks=element_blank(),axis.text = element_blank()) + # + #labs(x='',y='') +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"),axis.text =  element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  p3 <- plot_olm(doc_df2[Same_Inoculum==TRUE],'../Stat_Outputs/DOC_stats_H_Null.csv',med)  + 
    theme(axis.title=element_text(size=16)) + theme_bw() + #labs(x='',y='') +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  p4 <- plot_olm(doc_df2[Same_Inoculum==FALSE],'../Stat_Outputs/DOC_stats_I_Null.csv',med)  + 
    theme(axis.title=element_text(size=16)) + theme_bw() + #labs(x='',y='') +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
  p3_2 <- plot_olm(doc_df3[Same_Inoculum==TRUE],'../Stat_Outputs/DOC_stats_J_Null.csv',med2)  + 
    theme(axis.title=element_text(size=16)) + theme_bw() + #labs(x='',y='') +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  p4_2 <- plot_olm(doc_df3[Same_Inoculum==FALSE],'../Stat_Outputs/DOC_stats_K_Null.csv',med2)  + 
    theme(axis.title=element_text(size=16)) + theme_bw() + #labs(x='',y='') +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}
if(null==0){
  pMiddle = ggarrange(plot_olm2(doc_df2,null) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),labels='C')
  pMiddle_2 = ggarrange(plot_olm2_2(doc_df3,null) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),labels='B')
  pLeft <- ggarrange(p1 + annotation_custom(arrangeGrob(p1b),xmin=0.0,xmax=0.3,ymin=0.6,ymax=0.83),labels='B')
  pLeft_2 <- ggarrange(p1_2 + annotation_custom(arrangeGrob(p1b_2),xmin=0.0,xmax=0.3,ymin=0.6,ymax=0.83),labels='A')
}else{
  pMiddle = ggarrange(plot_olm2(doc_df2,null) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),labels='B')
  pMiddle_2 = ggarrange(plot_olm2_2(doc_df3,null) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()),labels='A')
  pLeft <- ggarrange(p1 + annotation_custom(arrangeGrob(p1b),xmin=0.0,xmax=0.3,ymin=0.6,ymax=0.83),labels='A')
  pLeft_2 <- ggarrange(p1_2 + annotation_custom(arrangeGrob(p1b_2),xmin=0.0,xmax=0.3,ymin=0.6,ymax=0.83),labels='A')
}

t1 = text_grob(expression('' == 'Inoculum'),face='bold',size=12)
t2 = text_grob(expression('' != 'Inoculum'),face='bold',size=12)
b = text_grob(expression('Overlap'),face='bold',size=16)
r = text_grob(expression('Dissimilarity'),face='bold',size=16,rot=90)
if(null==0){
  pRight = ggarrange(grid.arrange(p3,top=t1),grid.arrange(p4,top=t2),ncol=1,labels=c('D','E'))
  pRight_2 = ggarrange(grid.arrange(p3_2,top=t1),grid.arrange(p4_2,top=t2),ncol=1,labels=c('C','D'))
}else{pRight = ggarrange(grid.arrange(p3,top=t1),grid.arrange(p4,top=t2),ncol=1,labels=c('C','D'))
pRight_2 = ggarrange(grid.arrange(p3_2,top=t1),grid.arrange(p4_2,top=t2),ncol=1,labels=c('C','D'))
}
# 
f2 = grid.arrange(add_vp(pLeft),add_vp(pMiddle),add_vp(pRight),ncol=3,widths = c(2,1,1))
s2 = grid.arrange(add_vp(pLeft_2),add_vp(pMiddle_2),add_vp(pRight_2),ncol=3,widths = c(2,1,1))
if(null == 0){
  # save_plot('../Final Figures/Fig2.png',f2,base_height=13.5/2.8,base_width=27/2.8)
  save_plot('../Final Figures/Supp2.png',s2,base_height=13.5/2.8,base_width=27/2.8)}else{
  # save_plot('../Final Figures/Supp1.png',f2,base_height=13.5/2.8,base_width=27/2.8)
  save_plot('../Final Figures/Supp3.png',s2,base_height=13.5/2.8,base_width=27/2.8)
  }
# 
