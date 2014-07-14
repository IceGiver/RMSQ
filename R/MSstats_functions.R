#! /usr/bin/Rscript --vanilla --default-packages=utils
library(reshape2)
library(limma)
library(stats)
library(graphics)
library(sqldf)
library(ggplot2)
#library(MSstats.daily)
library(MSstats)
library(pheatmap)
library(RColorBrewer)
library(grDevices)
library(data.table)

source("~/Projects/HPCKrogan/Scripts/MSPipeline/src/Conversion/AnnotateWithUniprot_lib.R")
theme_set(theme_bw(base_size = 15, base_family="Helvetica"))

############
## FUNCTIONS

filterMaxqData = function(data){
  data_selected = data[grep("^CON__|^REV__",data$Proteins, invert=T),]
  return(data_selected)
}

explodeMaxQProteinGroups = function(data){
  return(data)
}

removeMaxQProteinGroups = function(data){
  data_selected = data[grep(";",data$Proteins, invert=T),]
  return(data_selected)
}

castMaxQToWide = function(d_long, aggregateFun=sum){
  data_w = dcast( Proteins + Sequence + Charge ~ Raw.file, data=d_long, value.var='Intensity', fun.aggregate=aggregateFun)
  return(data_w)
}

castMaxQToWidePTM = function(d_long, aggregateFun=sum){
  data_w = dcast( Proteins + Modified.sequence + Charge ~ Raw.file, data=d_long, value.var='Intensity', fun.aggregate=aggregateFun)
  setnames(data_w,2,'Sequence')
  return(data_w)
}

meltMaxQToLong = function(d_wide){
  data_l = melt(d_wide, id.vars=c('Proteins','Sequence','Charge'))
  setnames(data_l,old=4:5,new=c('Raw.file','Intensity'))
  return(data_l)
}

fillMissingMaxQEntries = function(data_w){
  mins = apply(data_w[,4:ncol(data_w)],2, function(x) min(x[x>0],na.rm=T))
  for(col in 4:ncol(data_w)){
    data_w[is.na(data_w[,col]),col]=0
    data_w[as.numeric(data_w[,col])==0,col]=mins[col-3]
  } 
  return(data_w)
}

mergeMaxQDataWithKeys = function(data_l, keys){
  unique_data = unique(data_l$Raw.file)
  unique_keys = unique(keys$RawFile)
  keys_not_found = setdiff(unique_keys, unique_data)
  data_not_found = setdiff(unique_data, unique_keys)
  cat(sprintf("keys found: %s \t keys not in data file:\n%s", length(unique_keys)-length(keys_not_found), paste(keys_not_found,collapse='\t')))
  cat(sprintf("data found: %s \t data not in keys file:\n%s", length(unique_data)-length(data_not_found), paste(data_not_found, collapse='\t')))
  
  ## select only required attributes from MQ format
  setnames(data_l, old='Raw.file', new='RawFile')
  data_l = merge(data_l, keys, by='RawFile')
  data_l
}

normalizePerCondition = function(d, NORMALIZATION_METHOD="scale"){
  unique_conditions = unique(d$Condition)
  d_tmp = c()
  
  for(u in unique_conditions){
    print(sprintf("normalizing\t%s",u))
    ss = data_w[d$Condition==u,]
    tmp = normalizeSingle(ss, NORMALIZATION_METHOD)
    d_tmp = rbind(d_tmp, tmp)
  }
  d_tmp
}

normalizeSingle = function(data_w, NORMALIZATION_METHOD="scale"){
  data_w_n = data.table(data_w[,1:3],normalizeBetweenArrays(as.matrix(data_w[,4:ncol(data_w)]), method=NORMALIZATION_METHOD))  
  data_w_n
}

dataToMSSFormat = function(d){
  tmp = data.frame(ProteinName=d$Proteins, PeptideSequence=d$Sequence, PrecursorCharge=NA, FragmentIon=NA, ProductCharge=d$Charge, IsotopeLabelType=d$IsotopeLabelType, Condition=d$Condition, BioReplicate=d$BioReplicate, Run=d$Run, Intensity=d$Intensity)
  tmp
}

###################################
## doesnt work so far with new code

normalizeToReference = function(d, keys, ref_protein=BAIT_REF, PDF=T, PRINT_DIR=""){
  
  ## merge background normalization info to data
  
  #d_l2 = merge(d, unique(keys[,c("BioReplicate","NormalizationGroup")]), by=c("BioReplicate"))
  d_l2=d
  ## log transform data
  d_l2$Intensity = log2(d$Intensity)  
  
  ########## ALL PROTEINS BEFORE 
  ##############################
  
  ## look at peptide distribution of all proteins
  if(PDF) pdf(sprintf("%s/AllProteinDistribution_b.pdf",PRINT_DIR), width=10, height=7)
  p = ggplot(data=d_l2, aes(x=BioReplicate, y=Intensity))
  print(p + geom_boxplot() +
          theme(axis.text.x=element_text(angle=-90)))  
  if(PDF) dev.off()
  
  ########## REF PROTEINS BEFORE 
  ##############################
  
  ## compute the average value for background reference peptides over all samples (background+bait)
  ## normalize against background reference proteins
  ref_peptides = d_l2[grep(refprotein, d_l2$Proteins),]
  ref_peptides = data.frame(ref_peptides, protein_id=as.character(paste0(ref_peptides$Sequence,ref_peptides$Charge)), stringsAsFactors=F)
  ref_peptides_counts = aggregate(BioReplicate ~ protein_id,data=ref_peptides,FUN=function(x)length(unique(x)))
  unique_replicates = length(unique(d_l2$BioReplicate))
  complete_observations = ref_peptides_counts[ref_peptides_counts$BioReplicate==unique_replicates,'protein_id']
  ref_peptides_complete = ref_peptides[ref_peptides$protein_id %in% complete_observations,]
  
  ## look at peptide distribution of ref proteins
  if(PDF) pdf(sprintf("%s/RefProteins_b.pdf",PRINT_DIR), width=10, height=7)
  p = ggplot(data=ref_peptides_complete, aes(x=BioReplicate, y=Intensity))
  print(p + geom_boxplot() +
          theme(axis.text.x=element_text(angle=-90)))  
  if(PDF) dev.off()
  
  ## look at their individual signal over replicates
  if(PDF) pdf(sprintf("%s/RefPeptides_b.pdf",PRINT_DIR), width=10, height=7)
  p = ggplot(data=ref_peptides_complete, aes(x=BioReplicate, y=Intensity, group=protein_id))
  print(p + geom_line(colour="grey") + 
          stat_smooth(aes(group=1), method = "lm",formula=y ~ poly(x, 6), colour="red", size=1) + 
          stat_summary(aes(group=1), geom="point", fun.y=median, shape=17, size=3, colour="darkred") + 
          ylab("Intensity") +
          ylim(15,30) +
          theme(axis.text.x=element_text(angle=-90)))
  if(PDF) dev.off()
  
  ref_peptides_avg_per_rep = aggregate(Intensity ~ BioReplicate, data=ref_peptides_complete, FUN=median)
  ref_peptides_avg_all = median(ref_peptides_avg_per_rep$Intensity)
  ref_peptides_avg_per_rep = data.frame(BioReplicate=ref_peptides_avg_per_rep$BioReplicate, correction=ref_peptides_avg_all-ref_peptides_avg_per_rep$Intensity)
  ref_peptides_complete = merge(ref_peptides_complete, ref_peptides_avg_per_rep, by='BioReplicate', all.x=T)
  ref_peptides_complete = data.frame(ref_peptides_complete, IntensityCorrected=ref_peptides_complete$Intensity+ref_peptides_complete$correction)
  
  ########## REF PROTEINS AFTER
  ##############################
  
  ## look at peptide distribution of ref proteins
  if(PDF) pdf(sprintf("%s/RefProteins_a.pdf",PRINT_DIR), width=10, height=7)
  p = ggplot(data=ref_peptides_complete, aes(x=BioReplicate, y=IntensityCorrected))
  print(p + geom_boxplot() +
          theme(axis.text.x=element_text(angle=-90)))  
  if(PDF) dev.off()
  
  ## look at their individual signal over replicates
  if(PDF) pdf(sprintf("%s/RefPeptides_a.pdf",PRINT_DIR), width=10, height=7)
  p = ggplot(data=ref_peptides_complete, aes(x=BioReplicate, y=IntensityCorrected, group=protein_id))
  print(p + geom_line(colour="grey") + 
          stat_smooth(aes(group=1), method = "lm",formula=y ~ poly(x, 6), colour="red", size=1) + 
          stat_summary(aes(group=1), geom="point", fun.y=median, shape=17, size=3, colour="darkred") + 
          ylab("Intensity") +
          ylim(15,30) +
          theme(axis.text.x=element_text(angle=-90)))
  if(PDF) dev.off()
  
  d_l2 = merge(d_l2, ref_peptides_avg_per_rep, by='BioReplicate', all.x=T)
  d_l2 = data.frame(d_l2, IntensityCorrected=d_l2$Intensity+d_l2$correction)
  
  ########## ALL PROTEINS BEFORE 
  ##############################
  
  ## look at peptide distribution of all proteins
  if(PDF) pdf(sprintf("%s/AllProteinDistribution_a.pdf",PRINT_DIR), width=10, height=7)
  p = ggplot(data=d_l2, aes(x=BioReplicate, y=IntensityCorrected))
  print(p + geom_boxplot() +
          theme(axis.text.x=element_text(angle=-90)))  
  if(PDF) dev.off()
  
  d_l2$Intensity = 2^d_l2$IntensityCorrected
  d_l2[,c("Run","Raw.file","Condition","BioReplicate","Proteins","Sequence","Charge","IsotopeLabelType","Intensity","NormalizationGroup")]
}
