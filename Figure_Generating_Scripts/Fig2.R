rm(list=ls())
library(data.table)
library(ggplot2)
library(ggpubr)
library(operators)

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

plot_high_overlap <- function(dat,med){
  odata = fread('../Data/Emergent_Comunity_Data_Equilibrium.csv')[Carbon_Source=='Glucose' & Transfer==12]
  cit_coms = unique(odata[ESV_ID == 'Citrobacter',]$Sample_ID)
  dat$Citrobacter <- 'One'
  dat[Var1 %!in% cit_coms & Var2 %!in% cit_coms,]$Citrobacter= 'Neither'
  dat[Var1 %in% cit_coms & Var2 %in% cit_coms,]$Citrobacter= 'Both'
  dat$Citrobacter = factor(dat$Citrobacter,levels=c('Both','One','Neither'))
  med =  floor(median(dat$Overlap)*100)/100
  dat = dat[Overlap>med,]
  pA <- ggplot() + geom_point(dat,mapping = aes(x=Overlap,y=Dissimilarity,col=Citrobacter)) +
    theme_pubr() + guides(fill=FALSE) + 
    # scale_x_continuous(breaks=c(floor(min(dat$Overlap)*1000)/1000,1),limits=c(floor(min(dat$Overlap)*1000)/1000,1)) + 
    scale_x_continuous(breaks=c(floor(med*100)/100,1),limits=c(floor(med*100)/100,1)) + 
    guides(alpha=FALSE,size=FALSE) +
    scale_y_continuous(breaks=c(0,0.8),limits=c(0,sqrt(log(2)))) +# stat_ellipse(type="norm",level=0.5) +
    theme(axis.title= element_text(size=16),axis.line = element_line(size=2),
          axis.text = element_text(size=14),legend.text = element_text(size=12)) + labs(col='Citrobacter In') +
    scale_color_manual(values=c('#0073C2FF','#8EE5EE','#EFC000FF')) +
    geom_point(dat[Var2 ==  'Glucose.2.6.12'],mapping = aes(x=Overlap,y=Dissimilarity),pch=21, fill=NA, size=6, colour="red")
  return(pA)
}

plot_comp <- function(dat,med,fn){
  odata = fread('../Data/Emergent_Comunity_Data_Equilibrium.csv')[Carbon_Source=='Glucose' & Transfer==12]
  cit_coms = unique(odata[ESV_ID == 'Citrobacter',]$Sample_ID)
  dat$Citrobacter <- 'One'
  dat[Var1 %!in% cit_coms & Var2 %!in% cit_coms,]$Citrobacter= 'Neither'
  dat[Var1 %in% cit_coms & Var2 %in% cit_coms,]$Citrobacter= 'Both'
  dat$Citrobacter = factor(dat$Citrobacter,levels=c('Both','One','Neither'))
  dat = dat[Overlap>med,]
  tests = fread(fn)
  tests = tests[Threshold ==med & !is.na(tests$t)]
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
  pvals$adj = paste('p < ', ceiling(1000/n)/1000)
  pvals[pvals$p !=0,]$adj = as.character(floor(pvals[pvals$p !=0,]$p*1000)/1000)
  pvals$group1 = c('Both','Both','One')
  pvals$group2 = c('One','Neither','Neither')
  pB <- ggboxplot(dat,x='Citrobacter',y='Dissimilarity',col='Citrobacter',palette = "jco",
                  add = "jitter",legend='right') +
    # stat_compare_means(comparisons = mycomparisons,method='t.test',size=2) + scale_y_continuous(breaks=c(0,0.8)) +
    guides(col=FALSE) +   labs(x = '') +
    stat_pvalue_manual(pvals,label='adj',y.position = c(0.8,0.9,0.85),size=3) +
    theme(legend.position = '' ,axis.title= element_text(size=16),axis.line = element_line(size=2),
          axis.text = element_text(size=14),legend.text = element_text(size=12)) 
  return(pB)
}

plot_one_pair<- function(){
  doc_df = fread('../Data/Dissimilarity_Overlap.csv')[Dissimilarity!=0 & Carbon_Source_1 == 'Glucose' & Carbon_Source_2 == 'Glucose' & Same_Inoculum==TRUE]
  odata = fread('../Data/Emergent_Comunity_Data.csv')[Carbon_Source == 'Glucose']
  doc_df = doc_df[(Var1 %in% gsub('.3','.12',unique(odata[Transfer==3]$Sample_ID))) & (Var2 %in% gsub('.3','.12',unique(odata[Transfer==3]$Sample_ID))),]
  odata$Label = paste('Inoculum ',odata$Inoculum,'Replicate ', odata$Replicate)
  
  cit_coms = unique(odata[ESV_ID == 'Citrobacter' & Transfer==12,]$Sample_ID)
  doc_df$Citrobacter <- 'One'
  doc_df[Var1 %!in% cit_coms & Var2 %!in% cit_coms,]$Citrobacter= 'Neither'
  doc_df[Var1 %in% cit_coms & Var2 %in% cit_coms,]$Citrobacter= 'Both'
  odata = odata[(Replicate %in% c(6,8) & Inoculum==2),]
  tax = unique(odata[Transfer>3 & Relative_Abundance>0.05]$ESV_ID)
  odata[-which(odata$ESV_ID %in% tax)]$ESV_ID = 'Other'
  odata = odata[,c('ESV_ID','Relative_Abundance','Label','Transfer','Family')]
  odata = odata[, lapply(.SD, sum, na.rm=TRUE), by=list(Label,ESV_ID,Transfer,Family), .SDcols=c("Relative_Abundance") ] 
  order_tax =   c('Enterobacteriaceae','Citrobacter','Raoultella','Aeromonas',
                  'Pseudomonas',"Pseudomonas.2",'Sphingobacterium','Other')
  odata$ESV_ID = factor(odata$ESV_ID,levels=order_tax)
  cm = c('blue','dodgerblue','lightblue','cadetblue',
         'firebrick','Red','darksalmon','orange','grey')
  p1<- Plot_TSeries(odata,cm)
  p2 <- frac_traj(odata)
  return(ggarrange(p1,p2,labels=c('E','F'),widths = c(2,2),common.legend=TRUE,legend='bottom'))
}


Plot_TSeries <- function(data,colours){
  p<- ggplot(data,aes(x=as.factor(Transfer),y=Relative_Abundance,fill=ESV_ID,width=0.85)) + 
    geom_bar(color='black',stat='identity') + 
    theme_classic() + 
    scale_fill_manual(drop=FALSE,values=colours) + #scale_fill_manual(values=colours_manual)
    theme(legend.position = "right")  + scale_x_discrete() + scale_y_continuous(breaks=c(0,0.5,1))  +
    labs(x='Transfer',y = 'Relative Abundance',fill='')   + 
    facet_wrap(~Label,nrow=2) + #guides(fill=FALSE)  +
    theme(axis.title = element_text(size=12),axis.text = element_text(size=10)
    ) + 
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )
  
  return(p)
}

frac_traj <- function(f_data){
  cit = f_data[grep('Citro',ESV_ID),]
  ent = f_data[grep('Enterob',Family),]
  ent = ent[, lapply(.SD, sum, na.rm=TRUE), by=list(Label,Transfer,Family), .SDcols=c("Relative_Abundance") ] 
  pseu1 = f_data[ESV_ID == 'Pseudomonas.2']
  pseu2 =  f_data[grep('Pseud',Family),]
  pseu2 = pseu2[, lapply(.SD, sum, na.rm=TRUE), by=list(Label,Transfer,Family), .SDcols=c("Relative_Abundance") ] 
  
  frac_a = c()
  frac_b = c()
  for(l in unique(f_data$Label)){
    for(t in 1:12){
      a = as.numeric(ent[Label==l & Transfer==t,'Relative_Abundance'])
      b =  as.numeric(cit[Label==l & Transfer==t,'Relative_Abundance'])
      b[is.na(b)] <-0
      a[is.na(a)] <-0
      frac_a = c(frac_a,b/(a))
      c = as.numeric(pseu1[Label==l & Transfer==t,'Relative_Abundance'])
      d =  as.numeric(pseu2[Label==l & Transfer==t,'Relative_Abundance'])
      c[is.na(c)] <-0
      d[is.na(d)] <-0
      frac_b = c(frac_b,c/d)
    }
  }
  traj_df = data.frame('Fent' = frac_a, 'Fpseud' = frac_b, Transfer = rep(1:12,2),Label = rep(unique(f_data$Label),each=12))
  traj_df$y_arrow = c(tail(traj_df$Fent,n=-1),NA)
  traj_df$x_arrow = c(tail(traj_df$Fpseud,n=-1),NA)
  traj_df$y_arrow[traj_df$Transfer == 12] <- NA
  traj_df$x_arrow[traj_df$Transfer  == 12] <- NA
  traj_df$Label = substr(as.character(traj_df$Label),13,nchar(as.character(traj_df$Label)))
  p <-   ggplot(traj_df, aes(x=Fpseud, y=Fent, col = Label)) + #geom_point() +
    geom_path()+ geom_segment(aes(xend=x_arrow, yend=y_arrow), arrow=arrow(length=unit(0.2,"cm")),size=0.5) + 
    theme_pubr()  + scale_colour_manual(values=c('Black','Red')) +  
    labs(col='',y=expression(frac('F'[C],'F'[Ent])),x=expression(frac('F'[P2],'F'[Pseud]))) +
    theme(legend.position = "right",
          legend.text = element_text(size=8),axis.line = element_line(size=1)) +
    scale_x_continuous(breaks=c(0,1)) + 
    scale_y_continuous(breaks=c(0,0.5)) + guides(col=FALSE) +
    theme(axis.title= element_text(size=12),axis.line = element_line(size=2),
          axis.text = element_text(size=14),legend.text = element_text(size=12)) 
  return(p)
}

doc_df = fread('../Data/Dissimilarity_Overlap.csv')
doc_df = doc_df[Dissimilarity!=0 & Same_Environment == TRUE & Carbon_Source_1 == 'Glucose',]
pA <- plot_high_overlap(doc_df[Same_Inoculum==TRUE],
                  med = median(doc_df[Same_Inoculum==TRUE]$Overlap))
pB <- plot_comp(doc_df[Same_Inoculum==TRUE],
                med = median(doc_df[Same_Inoculum==TRUE]$Overlap),
                fn = '../Stat_Outputs/TTest_CSource_Cit.csv')
pAB = ggarrange(pA,pB,labels = 'AUTO')
pCD = plot_univ(doc_df,'../Stat_Outputs/Curve_Fitting_Citrobacter.csv','../Stat_Outputs/Curve_Fitting_No_Citrobacter.csv')
pCD = ggarrange(pCD[[1]],pCD[[2]],ncol=2,labels=c('C','D'))
pEF = plot_one_pair()
ggsave('../Final_Figures/Fig2.png',grid.arrange(pAB,pCD,pEF,nrow=3,heights=c(6,6,6)),width=9,height=12)