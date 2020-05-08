rm(list=ls())
library(data.table)
library(reshape2)

set.seed(1)
cs = 'Glucose'
doc_df  = fread('../Data/Dissimilarity_Overlap.csv')
doc_df_b= fread('../Data/Dissimilarity_Overlap_Bootstrapped.csv')
doc_df = doc_df[Carbon_Source_1 == cs & Carbon_Source_2 == cs &Same_Inoculum ==TRUE & Dissimilarity!=0]
doc_df_b = doc_df_b[Carbon_Source_1 == cs & Carbon_Source_2 == cs &Same_Inoculum ==TRUE & Dissimilarity!=0]

med = median(doc_df$Overlap)
doc_df =doc_df[Overlap >med,]
doc_df = doc_df[,list(Overlap,Dissimilarity,Inoculum_1,Overlapping_Species)]

doc_df_b =doc_df_b[Overlap >med,]
doc_df_b = doc_df_b[,list(Overlap,Dissimilarity,Inoculum_1,Overlapping_Species,Run)]

odata = fread('../Data/Emergent_Comunity_Data_Equilibrium.csv')[Carbon_Source ==cs]
odata = aggregate(.~ESV_ID*Sample_ID,odata[,list(Relative_Abundance,ESV_ID,Sample_ID)],function(x) sum(x))


taxa_df =data.table(table(odata$ESV_ID))
taxa_df = taxa_df[rev(tail(order(taxa_df$N),10))]

for(j in taxa_df$V1){
  doc_df[,j := sapply(strsplit(doc_df$Overlapping_Species,'_'),function(x) j %in% x)]
  doc_df_b[,j := sapply(strsplit(doc_df_b$Overlapping_Species,'_'),function(x) j %in% x)]
  colnames(doc_df)[ncol(doc_df)] = j
  colnames(doc_df_b)[ncol(doc_df_b)] = j
}

# doc_df = doc_df[Citrobacter == FALSE,]
doc_df = doc_df[,!"Overlapping_Species"]
doc_df_b = doc_df_b[,!"Overlapping_Species"]

doc_df$Inoculum = as.factor(doc_df$Inoculum_1)
doc_df_b$Inoculum = as.factor(doc_df_b$Inoculum_1)

doc_df = doc_df[,!c('Inoculum_1','Overlap','Inoculum')]
doc_df_b = doc_df_b[,!c('Inoculum_1','Overlap','Inoculum')]

estimate = c()
pvalue =c()
pvalue2 =c()
statistic =c()
for(i in taxa_df$V1){
  t = t.test(doc_df[,Dissimilarity]~!as.vector(t(doc_df[,i,with=FALSE])),alternative='greater')
  statistic = c(statistic,t$statistic)
  estimate = c(estimate,t$estimate[1] - t$estimate[2])
  pvalue = c(pvalue,t$p.value[[1]])
  t_bootstrap =c()
  for(k in unique(doc_df_b$Run)){
    t_bootstrap = c(t_bootstrap,
                    tryCatch(t.test(doc_df_b[Run==k,Dissimilarity]~!as.vector(t(doc_df_b[Run==k,i,with=FALSE])),alternative='greater')$statistic,
                    error=function(e){return(NA)}))
  }
  t_bootstrap = t_bootstrap[!is.na(t_bootstrap)]
  pvalue2 = c(pvalue2,sum(t_bootstrap<0)/length(t_bootstrap))
}

taxa_df$'Difference in Mean' = estimate
taxa_df$'T-statistic' = statistic
taxa_df$'P.value' = as.character(pvalue2)
# taxa_df[taxa_df$P.value>1]$P.value = 1
taxa_df[taxa_df$P.value==0]$P.value = '<0.002'
colnames(taxa_df)[1] = 'ESV'
if(cs=='Glucose'){
  pdf("../Final_Figures/Table1.pdf", height=11, width=8.5)
  grid.table(taxa_df[rev(order(taxa_df$`T-statistic`)),], rows = NULL)
  dev.off()
}

if(cs=='Citrate'){
  pdf("../Final_Figures/TableS1.pdf", height=11, width=8.5)
  grid.table(taxa_df[rev(order(taxa_df$`T-statistic`)),], rows = NULL)
  dev.off()
}

if(cs=='Leucine'){
  pdf("../Final_Figures/TableS2.pdf", height=11, width=8.5)
  grid.table(taxa_df[rev(order(taxa_df$`T-statistic`)),], rows = NULL)
  dev.off()
}