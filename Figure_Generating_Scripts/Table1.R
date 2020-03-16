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
library(MASS)
library(stargazer)
set.seed(1)
ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 
'%!in%' <- function(x,y)!('%in%'(x,y))

OSR_calc <- function(x){
  x = as.character(x)
  return(length(strsplit(x,'_')[[1]]))
}

null = 0
if(null == 0){
  doc_df  = fread('../Data/Dissimilarity_Overlap.csv')
}else{
  doc_df = fread('../Data/Dissimilarity_Overlap_Null.csv')
}
cs = 'Glucose'
odata = fread('../Data/Emergent_Comunity_Data_Equilibrium.csv')[Carbon_Source ==cs]
doc_df = doc_df[Carbon_Source_1 == cs & Carbon_Source_2 == cs &Same_Inoculum ==TRUE]
doc_df = doc_df[Dissimilarity!=0,]
med = median(doc_df$Overlap)
print(med)
# odata = odata[Family=='Enterobacteriaceae']
odata = aggregate(.~ESV_ID*Sample_ID,odata[,list(Relative_Abundance,ESV_ID,Sample_ID)],function(x) sum(x))
taxa_df =data.table(table(odata$ESV_ID))
taxa_df = taxa_df[rev(tail(order(taxa_df$N),10))]
doc_df =doc_df[Overlap >med,]
doc_df = doc_df[,list(Overlap,Dissimilarity,Inoculum_1,Overlapping_Species)]
osl = strsplit(doc_df$Overlapping_Species,'_')
for(j in taxa_df$V1){
  doc_df[,j := sapply(osl,function(x) j %in% x)]
  colnames(doc_df)[ncol(doc_df)] = j
}
# doc_df = doc_df[Citrobacter == FALSE,]
doc_df = doc_df[,!"Overlapping_Species"]
doc_df$Inoculum = as.factor(doc_df$Inoculum_1)
doc_df = doc_df[,!c('Inoculum_1','Overlap','Inoculum')]
estimate = c()
pvalue =c()
statistic =c()
# taxa_df = taxa_df[taxa_df$V1 %!in% c('Citrobacter','Yersinia'),]
for(i in taxa_df$V1){
  t = t.test(doc_df[,Dissimilarity]~!as.vector(t(doc_df[,i,with=FALSE])),alternative='greater')
  statistic = c(statistic,t$statistic)
  estimate = c(estimate,t$estimate[1] - t$estimate[2])
  pvalue = c (pvalue,t$p.value[[1]])
}

taxa_df$'Difference in Mean' = estimate
taxa_df$'T-statistic' = statistic
taxa_df$'P.value' = pvalue *nrow(taxa_df)
taxa_df[taxa_df$P.value>1,]$P.value <- 1
colnames(taxa_df)[1] = 'ESV'
if(cs=='Glucose'){
  pdf("../Final Figures/Table1.pdf", height=11, width=8.5)
  grid.table(taxa_df[rev(order(taxa_df$`T-statistic`)),], rows = NULL)
  dev.off()
}

if(cs=='Citrate'){
  pdf("../Final Figures/TableS1.pdf", height=11, width=8.5)
  grid.table(taxa_df[rev(order(taxa_df$`T-statistic`)),], rows = NULL)
  dev.off()
}

if(cs=='Leucine'){
  pdf("../Final Figures/TableS2.pdf", height=11, width=8.5)
  grid.table(taxa_df[rev(order(taxa_df$`T-statistic`)),], rows = NULL)
  dev.off()
}