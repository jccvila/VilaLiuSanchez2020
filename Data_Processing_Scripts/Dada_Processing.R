rm(list=ls())
library(data.table)
library(readr)
set.seed(1)
#Script 1 for Extract data  from DADA2 output and rarefies to constant read depth. 
# Stores into a melted data.frame with standardized columns for subsequent analysis.

rarefy <- function(dat,n =min(colSums(dat))){
  # normalize to sample count
  for(i in 1:ncol(dat)){dat
    old_column = dat[,i,with=FALSE][[1]]
    elements = factor(row.names(dat),levels=row.names(dat))
    sample = table(sample(rep(elements,old_column),n,replace =FALSE))
    dat[,(i) := as.numeric(sample)]
  }
  return(dat)
}

assign_taxonomy_id <- function(dat){
  #Asigns ID to the higest available taxonomic level
  dat$ESV_ID = dat$Genus
  dat$ESV_ID[which(is.na(dat$ESV_ID))] = dat$Family[which(is.na(dat$ESV_ID))]
  dat$ESV_ID[which(is.na(dat$ESV_ID))] = dat$Order[which(is.na(dat$ESV_ID))]
  dat$ESV_ID[which(is.na(dat$ESV_ID))] = dat$Class[which(is.na(dat$ESV_ID))]
  dat$ESV_ID[which(is.na(dat$ESV_ID))] = dat$Phylum[which(is.na(dat$ESV_ID))]
  dat$ESV_ID[which(is.na(dat$ESV_ID))] = dat$Kingdom[which(is.na(dat$ESV_ID))]
  dat$ESV_ID = as.factor(make.unique(dat$ESV_ID))
  #Make NA use ESV_ID to avoid counting them as the same group
  dat[is.na(dat$Genus)]$Genus = as.character(dat[is.na(dat$Genus)]$ESV_ID)
  dat[is.na(dat$Family)]$Family = dat[is.na(dat$Family)]$Genus
  dat[is.na(dat$Order)]$Order = dat[is.na(dat$Order)]$Family
  dat[is.na(dat$Class)]$Class = dat[is.na(dat$Class)]$Order
  dat[is.na(dat$Phylum)]$Phylum = dat[is.na(dat$Phylum)]$Class
  
  return(dat)
}
set.seed(1) # To ensure reproducibility

#Load Data
aux = fread('../Data/metadata.csv')
Taxonomy_Data = fread('../Data/taxonomy.csv')
Abundance_Data = fread('../Data/otu_table.csv') #Data is actually at an ESV level

#Rarefy
Abundance_Data = as.matrix(Abundance_Data[,2:ncol(Abundance_Data)])
Abundance_Data = data.table(Abundance_Data)
Abundance_Data = rarefy(Abundance_Data) 

Taxonomy_Data = assign_taxonomy_id(Taxonomy_Data) #Assign Taxonomy ID
Abundance_Data = t(t(Abundance_Data)/colSums(Abundance_Data)) #Calculated Relative Abundance of ESVs



#Create Sample_ID from carbon source, community, replicate and trasnfer point
cref = aux$Carbon
cref[is.na(cref)] <- 'Original' #T0 inoculum
comref = parse_number(as.character(aux$Comm))
repref = as.numeric(aux$Rep)
tref = as.numeric(aux$Transfer)
colnames(Abundance_Data)= paste(cref,comref,repref,tref,sep='.')
Abundance_Data = data.table(Abundance_Data)
Abundance_Data$ESV_ID = Taxonomy_Data$ESV_ID

#Now convert matrix into a data.frame using melt
Abundance_Data = melt(Abundance_Data,id= c('ESV_ID'),variable.name ='Sample_ID',value.name = 'Relative_Abundance')
Abundance_Data = Abundance_Data[Relative_Abundance>0,]
#Extract caronb source, community,replicate number and transfer point from sample id

Abundance_Data$Carbon_Source =sapply(Abundance_Data$Sample_ID,function(x) strsplit(as.character(x),'[.]')[[1]][1])
Abundance_Data$Inoculum =sapply(Abundance_Data$Sample_ID,function(x) strsplit(as.character(x),'[.]')[[1]][2])
Abundance_Data$Replicate =sapply(Abundance_Data$Sample_ID,function(x) strsplit(as.character(x),'[.]')[[1]][3])
Abundance_Data$Transfer =sapply(Abundance_Data$Sample_ID,function(x) strsplit(as.character(x),'[.]')[[1]][4])
Abundance_Data$ESV_ID = droplevels(Abundance_Data$ESV_ID)

#Merge with taxonomy
merged_data = merge(Taxonomy_Data,Abundance_Data)

#Save Eveything
fwrite(merged_data,'../Data/Emergent_Comunity_Data.csv')
#Save Equilibrium
equilibrium = merged_data[Transfer==12]
fwrite(merged_data[Transfer==12],'../Data/Emergent_Comunity_Data_Equilibrium.csv')
