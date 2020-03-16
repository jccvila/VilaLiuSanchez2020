rm(list=ls())
library(data.table)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(ggpubr)
library(stringr)
library(stats)
library(grid)
set.seed(1)

odata = fread('../Data/Unnormalized_Emergent_Comunity_Data.csv')[Carbon_Source == 'Glucose']
'%!in%' <- function(x,y)!('%in%'(x,y))
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 
null = 0
if(null == 0){
  edf = fread('../Data/Dissimilarity_Overlap.csv')
  tdf = fread('../Data/Dissimilarity_Overlap_Tseries.csv')
}else{
  edf = fread('../Data/Dissimilarity_Overlap_Null.csv.csv')
  tdf = fread('../Data/Dissimilarity_Overlap_Tseries_Null.csv')
}

# doc_df = doc_df[Dissimilarity != 0 ,]
edf = edf[Same_Inoculum == TRUE & Same_Environment==TRUE & Carbon_Source_1 == 'Glucose']

tdf = fread('../Data/Dissimilarity_Overlap_Tseries.csv')
tdf = tdf[Same_Transfer==TRUE & Same_Inoculum==TRUE]

plot_Time<- function(data,c){
  xs = seq(0,1,length=21) + 0.025
  xs = xs[1:length(xs)-1]
  data= data[order(data$Transfer_1),]
  data$Label = factor(paste('T',data$Transfer_1),levels=unique(paste('T',data$Transfer_1)))
  p <- ggplot(data) + geom_point(aes(x=Overlap,y=Dissimilarity),col=c,alpha=1) +
    theme_bw() + guides(fill=FALSE) + 
    scale_x_continuous(breaks = c(0,1),limits=c(0,1)) +guides(col=FALSE) +
    scale_y_continuous(breaks=c(0,0.8),limits=c(0,sqrt(log(2)))) + facet_grid( ~ Label) + 
    theme(strip.text = element_text(size=10),axis.title = element_text(size=14),axis.text = element_text(size=12))
  return(p)
}



plot_Traj<- function(data,p){
  xs = seq(0,1,length=21) + 0.025
  xs = xs[1:length(xs)-1]
  data$pair = paste(substr(data$Var1,1,12),substr(data$Var2,1,12))
  data= data[order(data$Transfer_1),]
  data = data[order(data$pair)]
  data$x_arrow = c(tail(data$Overlap,n=-1),NA)
  data$y_arrow = c(tail(data$Dissimilarity,n=-1),NA)
  data$y_arrow[data$Transfer_1 == 12] <- NA
  data$x_arrow[data$Transfer_1 == 12] <- NA
  
  p <-   ggplot(data[pair==p,], aes(x=Overlap, y=Dissimilarity,col=pair)) + #geom_point() +
    geom_path()+ geom_segment(aes(xend=x_arrow, yend=y_arrow), arrow=arrow(length=unit(0.2,"cm")),size=0.5) + 
    theme_pubr()  + scale_colour_manual(values= c('Red')) + guides(col=FALSE) +
    labs(x='O',y='D') +
    theme(axis.line = element_line(size=1)) +
    scale_x_continuous(breaks=c(0,1),limits=c(0,1)) + 
    scale_y_continuous(breaks=c(0,1),limits=c(0,1)) + 
    theme(axis.title.y = element_text(size=12),axis.title.x=element_text(size=12))
  return(p)
}

Plot_TSeries <- function(data,colours){
  p<- ggplot(data,aes(x=as.factor(Transfer),y=Relative_Abundance,fill=ESV_ID,width=0.85)) + 
    geom_bar(color='black',stat='identity') + 
    theme_classic() + 
    scale_fill_manual(drop=FALSE,values=colours) + #scale_fill_manual(values=colours_manual)
    theme(legend.position = "right")  + scale_x_discrete() + scale_y_continuous(breaks=c(0,0.5,1))  +
    labs(x='Transfer',y = 'Relative Abundance',fill='')   + 
    facet_wrap(~Label,nrow=2,ncol=4) + #guides(fill=FALSE)  +
    theme(axis.title = element_text(size=12),axis.text = element_text(size=10)
    )
  
  return(p)
}



rep_no <-function(x){
  x = strsplit(x,'[.]')[[1]]
  return(x[3])
}

inoc_no <-function(x){
  x = strsplit(x,'[.]')[[1]]
  return(x[2])
}
frac_traj <- function(f_data){
  cit = f_data[grep('Citro',ESV_ID),]
  ent = f_data[grep('Entero',ESV_ID),]
  pseu1 = f_data[ESV_ID == 'Pseudomonas']
  pseu2 =  f_data[ESV_ID == 'Pseudomonas.2']
  frac_a = c()
  frac_b = c()
  for(l in unique(f_data$Label)){
    for(t in 1:12){
      a = as.numeric(ent[Label==l & Transfer==t,'Relative_Abundance'])
      b =  as.numeric(cit[Label==l & Transfer==t,'Relative_Abundance'])
      b[is.na(b)] <-0
      a[is.na(a)] <-0
      frac_a = c(frac_a,b/(a+b))
      c = as.numeric(pseu1[Label==l & Transfer==t,'Relative_Abundance'])
      d =  as.numeric(pseu2[Label==l & Transfer==t,'Relative_Abundance'])
      c[is.na(c)] <-0
      d[is.na(d)] <-0
      frac_b = c(frac_b,d/(c+d))
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
    labs(col='',y=expression(frac('F'[C],'F'[E] + 'F'[C])),x=expression(frac('F'[P],'F'[P] + 'F'[P2]))) +
    theme(legend.position = "right",
          legend.text = element_text(size=8),axis.line = element_line(size=1)) +
    scale_x_continuous(breaks=c(0,1)) + 
    scale_y_continuous(breaks=c(0,1)) + guides(col=FALSE) +
    theme(axis.title.y = element_text(size=12),axis.title.x=element_text(size=12))
  return(p)
}


plot_one_pair<- function(pair = 1){
  doc_df = fread('../Data/Dissimilarity_Overlap.csv')[Dissimilarity!=0 & Carbon_Source_1 == 'Glucose' & Carbon_Source_2 == 'Glucose' & Same_Inoculum==TRUE]
  odata = fread('../Data/Unnormalized_Emergent_Comunity_Data.csv')[Carbon_Source == 'Glucose']
  doc_df = doc_df[(Var1 %in% gsub('.3','.12',unique(odata[Transfer==3]$Sample_ID))) & (Var2 %in% gsub('.3','.12',unique(odata[Transfer==3]$Sample_ID))),]
  odata$Label = paste('Inoculum ',odata$Inoculum,'Replicate ', odata$Replicate)
  
  cit_coms = unique(odata[ESV_ID == 'Citrobacter' & Transfer==12,]$Sample_ID)
  doc_df$Citrobacter <- 'One'
  doc_df[Var1 %!in% cit_coms & Var2 %!in% cit_coms,]$Citrobacter= 'Neither'
  doc_df[Var1 %in% cit_coms & Var2 %in% cit_coms,]$Citrobacter= 'Both'
  pA <- ggplot(doc_df,aes(x=Overlap,y=Dissimilarity,col=Citrobacter)) + geom_point(size=3) +
    theme_pubr() + guides(fill=FALSE) + 
    # scale_x_continuous(breaks=c(floor(min(dat$Overlap)*1000)/1000,1),limits=c(floor(min(dat$Overlap)*1000)/1000,1)) + 
    scale_x_continuous(breaks=c(floor(med*1000)/1000,1),limits=c(floor(med*1000)/1000,1)) + 
    guides(alpha=FALSE,size=FALSE) +
    scale_y_continuous(breaks=c(0,0.8),limits=c(0,sqrt(log(2)))) +# stat_ellipse(type="norm",level=0.5) +
    theme(axis.title= element_text(size=16),axis.line = element_line(size=2),
          axis.text = element_text(size=14),legend.text = element_text(size=12)) + labs(col='Citrobacter In') +
    scale_color_manual(values=c('#0073C2FF','#8EE5EE','#EFC000FF'))
  
  
  
  
  
  odata = odata[(Replicate %in% c(6,8) & Inoculum==2),]
  tax = unique(odata[Transfer>3 & Relative_Abundance>0.05]$ESV_ID)
  odata[-which(odata$ESV_ID %in% tax)]$ESV_ID = 'Other'
  odata = odata[,c('ESV_ID','Relative_Abundance','Label','Transfer')]
  odata = odata[, lapply(.SD, sum, na.rm=TRUE), by=list(Label,ESV_ID,Transfer), .SDcols=c("Relative_Abundance") ] 
  order_tax =   c('Enterobacteriaceae','Citrobacter','Raoultella','Aeromonas',
                  'Pseudomonas',"Pseudomonas.2",'Sphingobacterium','Other')
  odata$ESV_ID = factor(odata$ESV_ID,levels=order_tax)
  cm = c('blue','dodgerblue','lightblue','cadetblue',
         'firebrick','Red','darksalmon','orange','grey')
  p1<- Plot_TSeries(odata,cm)
  p2 <- frac_traj(odata)
  tdf = fread('../Data/Dissimilarity_Overlap_Tseries.csv')
  tdf = tdf[Same_Transfer==TRUE & Same_Inoculum==TRUE & Inoculum_1==2]
  p3 <- plot_Traj(tdf,p="Glucose.2.8. Glucose.2.6.")
  ggsave('../Final Figures/Fig5.png',grid.arrange(ggarrange(p1,labels='A'),ggarrange(p2,p3,labels=c('B','C'),align='hv')),height=6,width=9)
  
}


plot_all_pairs <- function(){
  doc_df = fread('../Data/Dissimilarity_Overlap.csv')[Dissimilarity!=0 & Carbon_Source_1 == 'Glucose' & Carbon_Source_2 == 'Glucose' & Same_Inoculum==TRUE]
  odata = fread('../Data/Unnormalized_Emergent_Comunity_Data.csv')[Carbon_Source == 'Glucose']
  doc_df = doc_df[(Var1 %in% gsub('.3','.12',unique(odata[Transfer==3]$Sample_ID))) & (Var2 %in% gsub('.3','.12',unique(odata[Transfer==3]$Sample_ID))),]
  odata$Label = paste('Inoculum ',odata$Inoculum,'Replicate ', odata$Replicate)
  
  s= unique(gsub('.3','.12',odata[Transfer==3,]$Sample_ID))
  doc_df = doc_df[Var1 %in% s & Var2 %in% s,]
  
  
  odata = odata[(Replicate %in% c(6,8) & Inoculum==2),]
  tax = unique(odata[Transfer>3 & Relative_Abundance>0.05]$ESV_ID)
  odata[-which(odata$ESV_ID %in% tax)]$ESV_ID = 'Other'
  odata = odata[,c('ESV_ID','Relative_Abundance','Label','Transfer')]
  odata = odata[, lapply(.SD, sum, na.rm=TRUE), by=list(Label,ESV_ID,Transfer), .SDcols=c("Relative_Abundance") ] 
  order_tax =   c('Enterobacteriaceae','Citrobacter','Raoultella','Aeromonas',
                  'Pseudomonas',"Pseudomonas.2",'Sphingobacterium','Other')
  odata$ESV_ID = factor(odata$ESV_ID,levels=order_tax)
  cm = c('blue','dodgerblue','lightblue','cadetblue',
         'firebrick','Red','darksalmon','orange','grey')
  p1<- Plot_TSeries(odata,cm)
  p2 <- frac_traj(odata)
  tdf = fread('../Data/Dissimilarity_Overlap_Tseries.csv')
  tdf = tdf[Same_Transfer==TRUE & Same_Inoculum==TRUE & Inoculum_1==2]
  p3 <- plot_Traj(tdf,p="Glucose.2.8. Glucose.2.6.")
  ggsave('../Final Figures/Fig5.png',grid.arrange(ggarrange(p1,labels='A'),ggarrange(p2,p3,labels=c('B','C'),align='hv')),height=6,width=9)
  
}
add_vp <- function(x){return(grid.arrange(x,vp=viewport(width=0.8, height=0.9)))}

plot_one_pair(pair = 1)
