rm(list=ls())
library(data.table)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(ggExtra)
library(ggpubr)
library(stringr)
library(stats)
set.seed(1)

overlap <- function(i,j,data){
  #Calculate overlap between two relative abundance vectors
  x = data[,i]
  y = data[,j]
  return(sum(x[which(x>0 & y>0)] + y[which(x>0 &y>0)])/2)
}

disimilarity <- function(i,j,data) {
  #Calculate dissimilarity of overlapping species between two relative abundance vectors
  x = data[,i]
  y = data[,j]
  x_norm <- x[which(x>0 & y>0)]
  y_norm = y[which(x>0 & y>0)]
  x_norm <- x_norm/sum(x_norm)
  y_norm <- y_norm/sum(y_norm)
  m = (x_norm + y_norm)/2
  kld_x = sum(x_norm*log(x_norm/m))/2
  kld_y = sum(y_norm*log(y_norm/m))/2
  r = sqrt(kld_x+kld_y)
  if(is.nan(r)){r<-0}
  return(r)
}

Overlapping_Species <- function(i,j,data){
  #Determines the set of overlapping species between two samples
  i_data = data[Sample_ID == i,]
  j_data = data[Sample_ID == j,]

  return(paste(sort(intersect(i_data$ESV_ID,j_data$ESV_ID)),collapse='_'))
}

Overlapping_Genus <- function(i,j,data){
  #Determines the set of overlapping genuses between two samples
  i_data = data[Sample_ID == i,]
  j_data = data[Sample_ID == j,]
  return(paste(sort(intersect(i_data$Genus,j_data$Genus)),collapse='_'))
}

calc_overlap <- function(dat){
  #determines overlapping taxa between two communities
  ref = fread('../Data/Emergent_Comunity_Data.csv')
  dat[,Overlapping_Species := Overlapping_Species(Var1,Var2,ref),by=seq_len(nrow(dat))]
  dat[,Overlapping_Genus := Overlapping_Genus(Var1,Var2,ref),by=seq_len(nrow(dat))]
  
  return(dat)
}

rshuffle <- function(mat){
  #Randomizes the compositional data for the control analysis of DOC.  This corresponds to Null Model 1 in Bashan et al 2016
  mat = t(as.matrix(mat))
  for(i in 1:ncol(mat)){
    mat[mat[,i]>0,i] <- sample(mat[mat[,i]>0,i]) #Sample without replacement to shuffle
  }
  mat = mat/rowSums(mat)
  return(t(mat))
}


CS_text <- function(x){
  #Extract carbon source name from sample id
  return(strsplit(as.character(x),'[.]')[[1]][1])
}

In_text <- function(x){
  #EXtract inoculum from sample id
  return(strsplit(as.character(x),'[.]')[[1]][2])
}

doc_analysis <- function(dat,nullmethod){
  #Performs DOC analysis. Returns a list of sample ID, the overlap matrix and the dissimilarity matrix.
  s_order = unique(as.character(dat$Sample_ID))
  
  #convert metled data.frame back to realtive abundance matrix
  mat_form = dcast(data = dat,formula =ESV_ID~as.character(Sample_ID),fun.aggregate = sum,value.var = "Relative_Abundance")
  rownames(mat_form) = mat_form[,1]
  mat_form = as.matrix(mat_form[,-1])
  mat_form = mat_form[,match(s_order,colnames(mat_form))]
  
  if(nullmethod ==1){
    mat_form = rshuffle(mat_form)
  }
  n <-ncol(mat_form)
  
  #calculate overlap and dissimilarity between sample pairs
  ove <- Vectorize(overlap, vectorize.args=list("i","j"))
  dis<- Vectorize(disimilarity, vectorize.args=list("i","j"))
  ove_mat <- outer(1:n,1:n,ove,data=mat_form)
  dis_mat <- outer(1:n,1:n,dis,data=mat_form)
  
  sname = colnames(mat_form) #Sample names and assign as row and column names.
  rownames(ove_mat) =sname
  rownames(dis_mat) =sname
  colnames(ove_mat) =sname
  colnames(dis_mat) =sname
  ove_mat[upper.tri(ove_mat,diag=TRUE)] = NA
  dis_mat[upper.tri(dis_mat,diag=TRUE)] = NA
  return(list(sname,ove_mat,dis_mat))
}

format_doc_df <- function(sname,ove_mat,dis_mat){ 
  #Take output of DOC analysis and formats it into melted data.frame for subseuqent analysis.
  dis_df = melt(dis_mat,na.rm=F,value.name='Dissimilarity')
  doc_df = data.table(melt(ove_mat,na.rm=F,value.name='Overlap'))
  doc_df$Dissimilarity = dis_df$Dissimilarity 
  doc_df[,Carbon_Source_1 := sapply(Var1,CS_text)]
  doc_df[,Carbon_Source_2 := sapply(Var2,CS_text)]
  doc_df[,Inoculum_1 := sapply(Var1,In_text)]
  doc_df[,Inoculum_2 := sapply(Var2,In_text)]
  doc_df[,Same_Inoculum := (Inoculum_1 == Inoculum_2)]
  doc_df[,Same_Environment := (Carbon_Source_1 == Carbon_Source_2)]
  return(doc_df)
}

bootstrap_doc_df <- function(sname,ove_mat,dis_mat,runs){
  #Repeats DOC analysis on boot-strapped samples.
  bdoc_df  = format_doc_df(sname,ove_mat,dis_mat)
  bdoc_df$Run = 0
  for(j in 1:runs){
    print(j)
    set.seed(j)
    newsname = sample(sname,replace=TRUE) #Bootstrapp samples.
    new_ove_mat = ove_mat[newsname,newsname]
    new_dis_mat = dis_mat[newsname,newsname]
    tdf = format_doc_df(newsname,new_ove_mat,new_dis_mat)
    tdf$Run = j
    bdoc_df = rbind(bdoc_df,tdf)
  }
  return(bdoc_df[!is.na(Dissimilarity)])
}

set.seed(1)
#Load data and sort columns
dat = fread('../Data/Emergent_Comunity_Data_Equilibrium.csv')
dat = dat[order(dat$Replicate)]
dat = dat[order(dat$Inoculum)]
dat = dat[order(dat$Carbon_Source)]

#So calculate dissimilarity and overlap for main data and all bootstrapped samples.
null = 0
doc_results = doc_analysis(dat,null)
doc_df = format_doc_df(doc_results[[1]],doc_results[[2]],doc_results[[3]])
doc_df_b =bootstrap_doc_df(doc_results[[1]],doc_results[[2]],doc_results[[3]],500)
doc_df = doc_df[!is.na(doc_df$Overlap)]
doc_df_b = doc_df_b[!is.na(doc_df_b$Overlap)]
fwrite(doc_df,'../Data/Dissimilarity_Overlap.csv')
fwrite(doc_df_b,'../Data/Dissimilarity_Overlap_Bootstrapped.csv')
#So determine overlapping taxa for main data and all bootstrapped samples.
doc_df = calc_overlap(doc_df)
doc_df_b = calc_overlap(doc_df_b)
fwrite(doc_df,'../Data/Dissimilarity_Overlap.csv')
fwrite(doc_df_b,'../Data/Dissimilarity_Overlap_Bootstrapped.csv')

#repeat the exact same analysis for null samples
set.seed(1)
null = 1
doc_results = doc_analysis(dat,null)
doc_df = format_doc_df(doc_results[[1]],doc_results[[2]],doc_results[[3]])
doc_df_b =bootstrap_doc_df(doc_results[[1]],doc_results[[2]],doc_results[[3]],500)
doc_df = doc_df[!is.na(doc_df$Overlap)]
doc_df_b = doc_df_b[!is.na(doc_df_b$Overlap)]
fwrite(doc_df,'../Data/Dissimilarity_Overlap_Null.csv')
fwrite(doc_df_b,'../Data/Dissimilarity_Overlap_Bootstrapped_Null.csv')
doc_df = calc_overlap(doc_df)
doc_df_b = calc_overlap(doc_df_b)
fwrite(doc_df,'../Data/Dissimilarity_Overlap_Null.csv')
fwrite(doc_df_b,'../Data/Dissimilarity_Overlap_Bootstrapped_Null.csv')

