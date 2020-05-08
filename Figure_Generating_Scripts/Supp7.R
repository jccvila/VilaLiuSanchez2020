rm(list=ls())
library(data.table)
library(ggplot2)
library(ggpubr)


plot_lowess<- function(dat,fn,c){
  sdata= fread(fn)
  med = median(dat$Overlap)
  m = loess(dat$Dissimilarity ~ dat$Overlap,span=0.2,family="symmetric",iterations=5,statistics='none')
  sdata$y = predict(m,sdata$xs)
  sdata = sdata[!is.na(sdata$y)]
  sdata[sdata$LCI_lowess <0]$LCI_lowess = 0
  sdata[sdata$UCI_lowess>sqrt(log(2))]$UCI_lowess = sqrt(log(2))
  p<- ggplot(dat) +
    geom_point(data=dat,aes(x= Overlap,y = Dissimilarity),size=1,alpha=0.05) +
    theme_pubr()  +
    geom_segment(x=med,xend=med,y=0,yend=sqrt(log(2)),linetype=2,col=c) +
    geom_ribbon(data=sdata,aes(x=xs,ymin=LCI_lowess,ymax=UCI_lowess),alpha=0.5,fill=c) + 
    geom_line(data=sdata,aes(x=xs,y=y),col=c,size=2) +
    scale_x_continuous(breaks=c(0,1),limits=c(0,1)) +
    scale_y_continuous(breaks=c(0,0.8),limits=c(0,0.85))  + theme(axis.line = element_line(size=1))
  return(p)
}


plot_olm<- function(dat,fn,med,c){
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
    geom_point(data=dat,aes(x= Overlap,y = Dissimilarity),size=1,alpha=0.5) +
    theme_pubr() + 
    geom_line(data=sdata,aes(x=xs,y=mean_lm),col=c,size=2) +
    geom_ribbon(data=sdata,aes(x=xs,ymin=LCI_lm,ymax=UCI_lm),alpha=0.5,fill=c) +
    scale_x_continuous(breaks=c(floor(med*100)/100,1),limits=c(floor(med*100)/100,1)) +
    scale_y_continuous(breaks=c(0,0.8),limits=c(0,0.9)) + theme(axis.line = element_line(size=1)) +
    annotate(geom='text',x = floor(min(dat$Overlap)*100)/100,y=0.88,label=paste('\t\t\t\t m =', slope,' ' , pvalue),size=2.5)
  return(p)
}

plot_univ <- function(dat,fn1,fn2){
  odata = fread('../Data/Emergent_Comunity_Data_Equilibrium.csv')[Carbon_Source=='Glucose' & Transfer==12]
  cit_coms = unique(odata[ESV_ID == 'Citrobacter',]$Sample_ID)
  cit_dat = dat[Var1%in% cit_coms & Var2 %in% cit_coms]
  no_cit_dat = dat[Var1%!in% cit_coms & Var2 %!in% cit_coms]
  med1 = median(cit_dat$Overlap)
  med2 = median(no_cit_dat$Overlap)
  p2 <- plot_lowess(cit_dat,fn1,'#0073C2FF') +
    theme(axis.title= element_text(size=16),axis.line = element_line(size=2),
          axis.text = element_text(size=14))
  p2_inset <- plot_olm(cit_dat,fn1,med1,'#0073C2FF') +  theme_bw()  + 
    theme(axis.title=element_blank(),axis.ticks=element_blank(),axis.text = element_blank()) + # + #labs(x='',y='') +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"),axis.text =  element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p2 <- p2 + annotation_custom(arrangeGrob(p2_inset),xmin=0.0,xmax=0.35,ymin=0.58,ymax=0.83) 
  
  p3 <- plot_lowess(no_cit_dat,fn2,'#EFC000FF') +
    theme(axis.title= element_text(size=16),axis.line = element_line(size=2),
          axis.text = element_text(size=14))
  p3_inset <- plot_olm(no_cit_dat,fn2,med2,'#EFC000FF')+  theme_bw()  + 
    theme(axis.title=element_blank(),axis.ticks=element_blank(),axis.text = element_blank()) + # + #labs(x='',y='') +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"),axis.text =  element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p3 <- p3 + annotation_custom(arrangeGrob(p3_inset),xmin=0.0,xmax=0.35,ymin=0.58,ymax=0.83)
  
  tl = text_grob('Citrobacter Communities',face='bold',size=12)
  tr = text_grob('No Citrobacter Communities',face='bold',size=12)
  
  return(list(arrangeGrob(p2,top=tl),arrangeGrob(p3,top=tr)))
}





doc_df = fread('../Data/Dissimilarity_Overlap_Null.csv')
doc_df = doc_df[Dissimilarity!=0 & Same_Environment == TRUE & Carbon_Source_1 == 'Glucose',]
pAB = plot_univ(doc_df,'../Stat_Outputs/Curve_Fitting_Citrobacter_Null.csv','../Stat_Outputs/Curve_Fitting_No_Citrobacter_Null.csv')
pAB = ggarrange(pAB[[1]],pAB[[2]],ncol=2,labels=c('A','B'))
ggsave('../Final_Figures/Supp7.png',pAB,width=9,height=4)