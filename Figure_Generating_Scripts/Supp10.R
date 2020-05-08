rm(list=ls())
library(data.table)
library(ggplot2)

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



doc_df = fread('../Data/Dissimilarity_Overlap.csv')
doc_df = doc_df[Dissimilarity>0,]


add_vp <- function(x){return(grid.arrange(x,vp=viewport(width=0.9, height=0.9)))}

pA <- plot_lowess(doc_df,'../Stat_Outputs/Curve_Fitting_All.csv',median(doc_df$Overlap))
pAInset <- plot_olm(doc_df,'../Stat_Outputs/Curve_Fitting_All.csv',median(doc_df$Overlap)) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),axis.text =  element_blank(),axis.title = element_blank())
  
pB <- ggarrange(add_vp(plot_olm(doc_df[Same_Inoculum==TRUE & Same_Environment == TRUE],
               '../Stat_Outputs/Curve_Fitting_All_SameEnvSameInoc.csv',median(doc_df$Overlap)) +
  theme(axis.text =  element_blank(),axis.title = element_blank())),labels='B')

pC <- ggarrange(add_vp(plot_olm(doc_df[Same_Inoculum==FALSE & Same_Environment == TRUE],
               '../Stat_Outputs/Curve_Fitting_All_SameEnvDiffInoc.csv',median(doc_df$Overlap)) +
  theme(axis.text =  element_blank(),axis.title = element_blank())),labels='C')


  
pD <- ggarrange(add_vp(plot_olm(doc_df[Same_Inoculum==TRUE & Same_Environment == FALSE],
               '../Stat_Outputs/Curve_Fitting_All_DiffEnvSameInoc.csv',median(doc_df$Overlap)) +
  theme(axis.text =  element_blank(),axis.title = element_blank())),labels='D')


pE <- ggarrange(add_vp(plot_olm(doc_df[Same_Inoculum==FALSE & Same_Environment == FALSE],
               '../Stat_Outputs/Curve_Fitting_All_DiffEnvDiffInoc.csv',median(doc_df$Overlap)) +
  theme(axis.text =  element_blank(),axis.title = element_blank())),labels='E')

t1 = text_grob(expression('' == 'Inoculum'),face='bold',size=16,rot=90)
t2 = text_grob(expression('' != 'Inoculum'),face='bold',size=16,rot=90)
tl = text_grob(expression('' == 'Environment'),face='bold',size=16)
tr = text_grob(expression('' != 'Environment'),face='bold',size=16)
b = text_grob(expression('Overlap'),face='bold',size=16)
r = text_grob(expression('Dissimilarity'),face='bold',size=16,rot=90)
pLeft <- ggarrange(pA + annotation_custom(arrangeGrob(pAInset),xmin=0.0,xmax=0.3,ymin=0.6,ymax=0.83),labels='A')
pRight <- grid.arrange(tl,tr,t1,t2,pB,pC,pD,pE,b,r,layout_matrix=rbind(c(NA,NA,NA,1,NA,NA,NA,NA,2,NA,NA,NA),
                                                                       c(NA,5,5,5,5,5,7,7,7,7,7,NA),
                                                                       c(NA,5,5,5,5,5,7,7,7,7,7,NA),
                                                                       c(3,5,5,5,5,5,7,7,7,7,7,NA),
                                                                       c(NA,5,5,5,5,5,7,7,7,7,7,NA),
                                                                       c(NA,5,5,5,5,5,7,7,7,7,7,10),
                                                                       c(NA,6,6,6,6,6,8,8,8,8,8,10),
                                                                       c(NA,6,6,6,6,6,8,8,8,8,8,NA),
                                                                       c(4,6,6,6,6,6,8,8,8,8,8,NA),
                                                                       c(NA,6,6,6,6,6,8,8,8,8,8,NA),
                                                                       c(NA,6,6,6,6,6,8,8,8,8,8,NA),
                                                                       c(NA,NA,NA,NA,NA,9,9,NA,NA,NA,NA)))
ggsave(filename = '../Final_Figures/Supp10.png',ggarrange(pLeft,pRight),height=5,width=10)