rm(list=ls())
library(data.table)
library(readr)
set.seed(1)



#Download Data
aux = fread('../Data/metadata.csv')
Taxonomy_Data = fread('../Data/taxonomy.csv')
OTU_Data = fread('../Data/otu_table.csv')
OTU_Data = as.matrix(OTU_Data[,2:ncol(OTU_Data)])
rownames(OTU_Data) = seq(1,nrow(OTU_Data))
colnames(OTU_Data) = seq(1,ncol(OTU_Data))
OTU_Data = data.table(OTU_Data)
sample_size = min(colSums(OTU_Data))

# normalize to sample count
for(i in 1:ncol(OTU_Data)){OTU_Data
  old_column = OTU_Data[,i,with=FALSE][[1]]
  elements = factor(row.names(OTU_Data),levels=row.names(OTU_Data))
  sample = table(sample(rep(elements,old_column),sample_size,replace =FALSE))
  OTU_Data[,(i) := as.numeric(sample)]
}

#Extract Taxonomic Data and labels Familys (I'll create an ID for each ESV based on the highest taxonomic level that  has been assigned
OTU_Data$ESV_ID = Taxonomy_Data$Genus
OTU_Data$ESV_ID[which(is.na(OTU_Data$ESV_ID))] = Taxonomy_Data$Family[which(is.na(OTU_Data$ESV_ID))]
OTU_Data$ESV_ID[which(is.na(OTU_Data$ESV_ID))] = Taxonomy_Data$Order[which(is.na(OTU_Data$ESV_ID))]
OTU_Data$ESV_ID[which(is.na(OTU_Data$ESV_ID))] = Taxonomy_Data$Class[which(is.na(OTU_Data$ESV_ID))]
OTU_Data$ESV_ID[which(is.na(OTU_Data$ESV_ID))] = Taxonomy_Data$Phylum[which(is.na(OTU_Data$ESV_ID))]
OTU_Data$ESV_ID[which(is.na(OTU_Data$ESV_ID))] = Taxonomy_Data$Kingdom[which(is.na(OTU_Data$ESV_ID))]
OTU_Data$ESV_ID = as.factor(make.unique(OTU_Data$ESV_ID))
Taxonomy_Data$ESV_ID = OTU_Data$ESV_ID


#Calculated Relative Abundance at a given taxonomic level (here ESV )
OTU_Data = OTU_Data[ ,lapply(.SD,sum), by = ESV_ID]
colnames(OTU_Data) = c('ESV_ID',as.character(aux$Carbon))
OTU_Data= OTU_Data[rowSums(OTU_Data[,2:ncol(OTU_Data)])>0.0,]
RAD = OTU_Data[,2:ncol(OTU_Data)]
RAD = data.table(t(t(as.matrix(RAD))/colSums(as.matrix(RAD))))

#Create ID 
cref = colnames(RAD)
comref = parse_number(as.character(aux$Comm))
repref = as.numeric(aux$Rep)
tref = as.numeric(aux$Transfer)
colnames(RAD)= paste(cref,comref,repref,tref,sep='.')
RAD$ESV_ID = OTU_Data$ESV_ID

#Now convert matrix into a data.frame using melt
plot_OTU_data = melt(RAD,id= c('ESV_ID'))
plot_OTU_data = plot_OTU_data[value>0,]
plot_OTU_data$Carbon_Source =sapply(plot_OTU_data$variable,function(x) strsplit(as.character(x),'[.]')[[1]][1])
plot_OTU_data$Com_No =sapply(plot_OTU_data$variable,function(x) strsplit(as.character(x),'[.]')[[1]][2])
plot_OTU_data$Replicate_No =sapply(plot_OTU_data$variable,function(x) strsplit(as.character(x),'[.]')[[1]][3])
plot_OTU_data$Transfer_No =sapply(plot_OTU_data$variable,function(x) strsplit(as.character(x),'[.]')[[1]][4])
plot_OTU_data$ESV_ID = droplevels(plot_OTU_data$ESV_ID)
colnames(plot_OTU_data) =c('ESV_ID','Sample_ID','Relative_Abundance','Carbon_Source','Inoculum','Replicate','Transfer')

#Nan are the original inoculum
plot_OTU_data[grep('NaN',plot_OTU_data$Sample_ID),]$Sample_ID <- paste('Original',sapply(plot_OTU_data[grep('NaN',plot_OTU_data$Sample_ID),]$Sample_ID,function(x) strsplit(as.character(x),'[.]')[[1]][2]),sep='.')
plot_OTU_data[grep('Original',plot_OTU_data$Sample_ID),]$Carbon_Source = 'NaN'

#Merge with taxonomy
merged_data = merge(Taxonomy_Data,plot_OTU_data)
merged_data[is.na(merged_data$Genus)]$Genus = as.character(merged_data[is.na(merged_data$Genus)]$ESV_ID)
merged_data[is.na(merged_data$Family)]$Family = merged_data[is.na(merged_data$Family)]$Genus
merged_data[is.na(merged_data$Order)]$Order = merged_data[is.na(merged_data$Order)]$Family
merged_data[]
#Save Eveything
fwrite(merged_data,'../Data/Emergent_Comunity_Data.csv')
#Save Equilibrium
equilibrium = merged_data[Transfer==12]
fwrite(merged_data[Transfer==12],'../Data/Emergent_Comunity_Data_Equilibrium.csv')
