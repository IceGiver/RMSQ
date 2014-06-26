#! /usr/bin/Rscript --vanilla --default-packages=utils
library(reshape2)
library(limma)
library(stats)
library(graphics)
library(sqldf)
library(ggplot2)
#library(MSstats.daily)
library(MSstats)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(grDevices)
source("~/Projects/HPCKrogan/Scripts/MSPipeline/src/Conversion/AnnotateWithUniprot_lib.R")

############
## FUNCTIONS

padMissing = function(data_mss_format){
  data_mss_format = data.frame(data_mss_format, stringsAsFactors=F)
  ## GET UNIQUE OBSERVATIONS ACROSS ALL CONDITIONS 
  unique_columns = c("ProteinName", "PeptideSequence", "PrecursorCharge","FragmentIon", "ProductCharge", "IsotopeLabelType")
  unique_observations = unique(data_mss_format[,unique_columns])
  ## GET UNIQUE CONDITIONS
  exp_columns = c("Condition","BioReplicate","Run")
  exp_set = unique(data_mss_format[, exp_columns])
  unique_observations = data.frame(id=1:nrow(unique_observations), unique_observations, stringsAsFactors=F)
  data_mss_indexed = merge(data_mss_format, unique_observations, by=unique_columns, stringsAsFactors=F)
  data_mss_indexed$Run = as.character(data_mss_indexed$Run)
  unique_ids = unique_observations$id
  unique_exps = unique(data_mss_format$BioReplicate)
  missing_data = c()
  
  for(r in unique_exps){
    run_sub_set = data_mss_indexed[data_mss_indexed$BioReplicate == r,]
    run_id = run_sub_set$Run[1] 
    present_ids = run_sub_set$id
    missing_ids = unique_ids[ ! unique_ids %in% present_ids]
    if(length(missing_ids)>0){
      missing_data = rbind(missing_data, cbind(run_id, missing_ids, NA)) 
    }
  }
  if(!is.null(missing_data)){
    colnames(missing_data) = c("Run", "id", "Intensity")
    total_data = rbind(data_mss_indexed[,c("Run", "id", "Intensity")], missing_data) 
    
    data_mss_padded = merge(exp_set, total_data, by=c("Run")) 
    data_mss_padded = merge(data_mss_padded, unique_observations, by=c("id")) 
    data_mss_padded = data_mss_padded[,colnames(data_mss_format)]
    #print(sqldf("select Run, count(*) from data_mss_padded group by Run"))
    data_mss_padded  
  }else{
    data_mss_format
  }
}

fillMissing = function(dw, rng_data){
  ## replace NA with minimum observed value 
  mins = apply(dw, 2, min, na.rm=T)
  tmp=dw
  for(col in rng_data){
    print(sprintf("filling %s with %s",colnames(dw)[col],mins[[col]]))
    tmp[is.na(tmp[,col]),col]=as.numeric(mins[[col]])
  }
  tmp
}

### READ DATA ####################################

maxqFilterAndJoin = function(data, keys){
  ## get rid of CONTAMINANTS
  data_selected = data[grep("CON__|REV__|;",data$Proteins, invert=T),]
  ## select only required attributes from MQ format
  data_selected = data_selected[, c("Proteins","Sequence","Charge","Raw.file","Intensity")]
  data_selected = data_selected[!is.na(data_selected$Intensity),]
  data_selected_keys = merge(data_selected, keys, by.x="Raw.file", by.y="RawFile")
  tmp = unique(data_selected_keys$Run)
  keys_not_found = paste(setdiff(keys$Run, tmp),sep=',')
  if(length(keys_not_found)==0){
    keys_not_found = 'none'
  }
  print(sprintf("keys not found in data file: %s", keys_not_found))
  data_selected_keys
}

maxqToWide = function(data){
  ## convert to wide format
  data_selected_wide = dcast(data=data, formula=Proteins + Sequence + Charge ~ Run + IsotopeLabelType, value.var="Intensity",fun.aggregate=median, na.rm=T)
  data_selected_wide
}

plotIntensityData = function(d){
  p = ggplot(data=d, aes(x=Run, y=log2(Intensity)))
  print(p + geom_boxplot() + facet_wrap(~ Condition, scales="free_x") + theme(axis.text.x=element_text(angle=45, hjust=1)))
}

theme_set(theme_bw(base_size = 15, base_family="Helvetica"))

normalizeToProtein = function(d, keys, bckgr_ref=BACKGROUND_REF, bait_ref=BAIT_REF, PDF=T, PRINT_DIR="Summary/20140403/"){
  
  ## merge background normalization info to data
  
  #d_l2 = merge(d, unique(keys[,c("BioReplicate","NormalizationGroup")]), by=c("BioReplicate"))
  d_l2=d
  ## log transform data
  d_l2$Intensity = log2(d$Intensity)  
  
  ## look at peptide distribution of all proteins
  if(PDF) pdf(sprintf("%s/AllProteinDistribution.pdf",PRINT_DIR), width=10, height=7)
  p = ggplot(data=d_l2, aes(x=BioReplicate, y=Intensity))
  print(p + geom_boxplot() +
    theme(axis.text.x=element_text(angle=-90)))  
  if(PDF) dev.off()
  
  ## get reference proteins to normalize against
  ref_proteins = paste(c(BACKGROUND_REF, BAIT_REF), collapse="|")
  ref_peptides = d_l2[grep(ref_proteins, d_l2$Proteins),]
  ref_peptides = data.frame(ref_peptides, reference=ifelse(grepl(BAIT_REF,ref_peptides$Proteins),"BAIT","BACKGROUND"))
  ref_peptides = data.frame(ref_peptides, protein_id=paste0(ref_peptides$Sequence,ref_peptides$Charge))
  
  ## look at peptide distribution of ref proteins
  if(PDF) pdf(sprintf("%s/RefProteinDistribution.pdf",PRINT_DIR), width=10, height=7)
  p = ggplot(data=ref_peptides, aes(x=BioReplicate, y=Intensity))
  print(p + geom_boxplot() +
    theme(axis.text.x=element_text(angle=-90)))  
  if(PDF) dev.off()
  
  ## look at their individual signal over replicates
  if(PDF) pdf(sprintf("%s/RefProteinPeptideSignal.pdf",PRINT_DIR), width=10, height=7)
  p = ggplot(data=ref_peptides, aes(x=BioReplicate, y=Intensity, group=protein_id))
  print(p + geom_line(colour="grey") + 
    stat_smooth(aes(group=1), method = "lm",formula=y ~ poly(x, 6), colour="red", size=1) + 
    stat_summary(aes(group=1), geom="point", fun.y=mean, shape=17, size=3, colour="darkred") + 
    ylab("Intensity") +
    ylim(15,30) +
    theme(axis.text.x=element_text(angle=-90)) + 
    facet_wrap( ~reference, scales="free_x"))
  if(PDF) dev.off()
  
  ## compute the average value for background reference peptides over all samples (background+bait)
  ref_peptides_bckg = ref_peptides[(ref_peptides$reference == 'BACKGROUND'), ] 
  ref_peptides_bckg_avg = aggregate(Intensity ~ Proteins+Sequence+Charge+reference, data=ref_peptides_bckg, FUN=AGGREGATE_FUN)
  colnames(ref_peptides_bckg_avg)[5] = "Intensity_avg"
  ref_peptides_bckg_adjm = merge(ref_peptides_bckg, ref_peptides_bckg_avg, by=c("Proteins", "Sequence", "Charge","reference"))
  ref_peptides_bckg_adjm = data.frame(ref_peptides_bckg_adjm, bckg_scale=ref_peptides_bckg_adjm$Intensity_avg/ref_peptides_bckg_adjm$Intensity)
  ref_peptides_bckg_adjm_avg = aggregate( bckg_scale ~ BioReplicate, data=ref_peptides_bckg_adjm, FUN=AGGREGATE_FUN)
  d_adj_bckg = merge(d_l2, ref_peptides_bckg_adjm_avg[,c("BioReplicate","bckg_scale")], by=c("BioReplicate"))
  d_adj_bckg = data.frame(d_adj_bckg, Intensity_bckg_adj=d_adj_bckg$Intensity*d_adj_bckg$bckg_scale)
  
  ## get the reference peptides another time
  ref_peptides = d_adj_bckg[grep(ref_proteins, d_adj_bckg$Proteins),]
  ref_peptides = data.frame(ref_peptides, reference=ifelse(grepl(BAIT_REF,ref_peptides$Proteins),"BAIT","BACKGROUND"))
  ref_peptides = data.frame(ref_peptides, protein_id=paste0(ref_peptides$Sequence,ref_peptides$Charge))
  
  ## look at their individual signal over replicates
  if(PDF) pdf(sprintf("%s/RefProteinPeptideSignal_bckg_normalized.pdf",PRINT_DIR), width=10, height=7)
  p = ggplot(data=ref_peptides, aes(x=BioReplicate, y=Intensity_bckg_adj, group=protein_id))
  print(p + geom_line(colour="grey") + 
    stat_smooth(aes(group=1), method = "lm",formula=y ~ poly(x, 6), colour="red", size=1) + 
    stat_summary(aes(group=1), geom="point", fun.y=mean, shape=17, size=3, colour="darkred") + 
    ylab("Intensity") +
    ylim(15,30) +
    theme(axis.text.x=element_text(angle=-90)) + 
    facet_wrap( ~reference, scales="free_x"))
  if(PDF) dev.off()
  
  ref_peptides_bait = ref_peptides[(ref_peptides$reference == 'BAIT' & ref_peptides$NormalizationGroup == 'BAIT'), ] 
  ref_peptides_bait_avg = aggregate(Intensity_bckg_adj ~ Proteins+Sequence+Charge+reference, data=ref_peptides_bait, FUN=AGGREGATE_FUN)
  colnames(ref_peptides_bait_avg)[5] = "Intensity_avg"
  ref_peptides_bait_adjm = merge(ref_peptides_bait, ref_peptides_bait_avg, by=c("Proteins", "Sequence", "Charge", "reference"))
  ref_peptides_bait_adjm = data.frame(ref_peptides_bait_adjm, bait_scale=ref_peptides_bait_adjm$Intensity_avg/ref_peptides_bait_adjm$Intensity_bckg_adj)
  ref_peptides_bait_adjm_avg = aggregate( bait_scale ~ BioReplicate, data=ref_peptides_bait_adjm, FUN=AGGREGATE_FUN)
  d_adj_bait = merge(d_adj_bckg, ref_peptides_bait_adjm_avg[,c("BioReplicate","bait_scale")], by=c("BioReplicate"), all.x=T)
  d_adj_bait = data.frame(d_adj_bait, Intensity_bait_adj=ifelse(d_adj_bait$NormalizationGroup =="BACKGROUND", d_adj_bait$Intensity_bckg_adj, d_adj_bait$Intensity_bckg_adj * d_adj_bait$bait_scale))
  
  d_adj = d_adj_bait
  
  ## get reference proteins to check normalization effect
  ref_peptides = d_adj[grep(ref_proteins, d_adj$Proteins),]
  ref_peptides = data.frame(ref_peptides, reference=ifelse(grepl(BAIT_REF,ref_peptides$Proteins),"BAIT","BACKGROUND"))
  ref_peptides = data.frame(ref_peptides, protein_id=paste0(ref_peptides$Sequence,ref_peptides$Charge))
  
  ## look at their individual signal over replicates after normalization
  if(PDF) pdf(sprintf("%s/RefProteinPeptideSignal_bck_bait_normalized.pdf",PRINT_DIR), width=10, height=7)
  p = ggplot(data=ref_peptides, aes(x=BioReplicate, y=Intensity_bait_adj, group=protein_id))
  print(p + geom_line(colour="grey") + 
    stat_smooth(aes(group=1), method = "lm",formula=y ~ poly(x, 6), colour="red", size=1) + 
    stat_summary(aes(group=1), geom="point", fun.y=mean, shape=17, size=3, colour="darkred") + 
    ylab("Intensity") +
    ylim(15,30) +
    theme(axis.text.x=element_text(angle=-90)) + 
    facet_wrap( ~reference, scales="free_x"))
  if(PDF) dev.off()
  
#     ## plot the scaling adjustement factors
#     if(PDF) pdf(sprintf("%sRefProteinsScaling.pdf",PRINT_DIR), width=10, height=7)
#     ggplot(ref_peptides_rep_avg, aes(y=1/scale, x=BioReplicate, fill=NormalizationGroup)) + 
#       geom_bar(stat="identity", position="dodge") +
#       theme(axis.text.x=element_text(angle=-90)) + 
#       ylab("Scale Adjustment factor") +
#       geom_hline(aes(yintercept=1))
#     if(PDF) dev.off()
#   
  ## look at peptide distribution of ref proteins
  if(PDF) pdf(sprintf("%s/RefProteinDistribution_normalized.pdf",PRINT_DIR), width=10, height=7)
  p = ggplot(data=ref_peptides, aes(x=BioReplicate, y=Intensity_bait_adj))
  print(p + geom_boxplot() +
    theme(axis.text.x=element_text(angle=-90)))
  if(PDF) dev.off()
  
  ## look at peptide distribution of all proteins
  if(PDF) pdf(sprintf("%s/AllProteinDistribution_normalized.pdf",PRINT_DIR), width=10, height=7)
  p = ggplot(data=d_adj, aes(x=BioReplicate, y=Intensity_bait_adj))
  print(p + geom_boxplot() +
    theme(axis.text.x=element_text(angle=-90)))  
  if(PDF) dev.off()

  d_adj$Intensity = 2^d_adj$Intensity_bait_adj
  d_adj[,c("Run","Raw.file","Condition","BioReplicate","Proteins","Sequence","Charge","IsotopeLabelType","Intensity")]
}

normalizeSingle = function(d, NORMALIZATION_METHOD="scale", FILL_MISSING=F, ref_name='HCVJF_NS5A'){
  unique_key_cols = unique(d[,c("Raw.file","Condition","BioReplicate","Run")])
  unique_val_cols = c("Proteins","Sequence","Charge")
  
  d_w = maxqToWide(d)
  if(NORMALIZATION_METHOD=="toref"){ ##custom
    d_w_n = normalizeToRef(d_w, bait=bait_name)
  }else{
    d_w_n = data.frame(d_w[,1:3],normalizeBetweenArrays(as.matrix(d_w[,4:ncol(d_w)]), method=NORMALIZATION_METHOD))  
  }
  if(FILL_MISSING){
    print("FILLING WITH MINIMUM VALUES")
    d_w_n = fillMissing(dw=d_w_n,rng_data=4:ncol(d_w_n))
  }
  write.table(d_w_n, eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
  d_w_n_m = melt(d_w_n, id.vars=c("Proteins","Sequence","Charge"), value.name="Intensity",variable.name='Run',na.rm=F)
  d_w_n_m_rf = data.frame(d_w_n_m[,unique_val_cols],IsotopeLabelType=gsub('([0-9,A-Z,\\.]+)_([A-Z,0-9]+)','\\2',d_w_n_m$Run),Run=gsub('([0-9,A-Z,\\.]+)_([A-Z,0-9]+)','\\1',d_w_n_m$Run),Intensity=d_w_n_m$Intensity)
  tmp = merge(unique_key_cols, d_w_n_m_rf, by="Run")
  tmp
}

normalizePerCondition = function(d, NORMALIZATION_METHOD="scale",  FILL_MISSING=F){
  unique_conditions = unique(d$Condition)
  d_tmp = c()
  
  for(u in unique_conditions){
    print(sprintf("normalizing\t%s",u))
    ss = d[d$Condition==u,]
    tmp = normalizeSingle(ss, NORMALIZATION_METHOD, FILL_MISSING)
    d_tmp = rbind(d_tmp, tmp)
  }
  d_tmp
}

sampleProteins = function(d, SAMPLE_SIZE=1000){
  unique_conditions = unique(d$Condition)
  unique_proteins = unique(d$Proteins)
  if(length(unique_proteins)<SAMPLE_SIZE) SAMPLE_SIZE=length(unique_proteins)
  protein_sample = sample(unique_proteins, SAMPLE_SIZE)
  
  d_tmp = c()
  for(u in unique_conditions){
    ss = d[d$Condition==u,]
    tmp = ss[ss$Proteins %in% protein_sample,]
    d_tmp = rbind(d_tmp, tmp)
  }
  d_tmp
}

dataToMSSFormat = function(d){
  tmp = data.frame(ProteinName=d$Proteins, PeptideSequence=d$Sequence, PrecursorCharge=NA, FragmentIon=NA, ProductCharge=d$Charge, IsotopeLabelType=d$IsotopeLabelType, Condition=d$Condition, BioReplicate=d$BioReplicate, Run=d$Run, Intensity=d$Intensity)
  tmp
}

uniprot_to_gene_replace = function(uniprot_ac){
  GENES = read.delim("~/Projects/HPCKrogan/Scripts/MSPipeline/files/uniprot_id_to_genes.txt", stringsAsFactors=F)
  merged = merge(data.frame(uniprot_ac=uniprot_ac, stringsAsFactors=F), GENES, all.x=T, by.x="uniprot_ac", by.y="uniprot_id")
  merged[is.na(merged$gene_name),]$gene_name = as.character(merged[is.na(merged$gene_name),'uniprot_ac'])
  merged
}

plotHeat = function(mss_F, out_file){
  heat_data = data.frame(mss_F, names=paste(mss_F$Protein,mss_F$uniprot_id,sep=' | '))
  heat_data = heat_data[,c('names','Label','log2FC')]
  heat_data_w = dcast(names ~ Label, data=heat_data, value.var='log2FC')
  #gene_names = uniprot_to_gene_replace(uniprot_ac=heat_data_w$Protein)
  rownames(heat_data_w) = heat_data_w$names
  heat_data_w = heat_data_w[,-1]
  heat_data_w[is.na(heat_data_w)]=0
  max_val = ceiling(max(heat_data_w))
  min_val = floor(min(heat_data_w))
  range_val = abs(min_val)+max_val
  bin_size=1
  signed_bins = (range_val/bin_size/2)-1
  colors_neg = rev(colorRampPalette(brewer.pal("Blues",n=5))(signed_bins))
  colors_pos = colorRampPalette(brewer.pal("Reds",n=5))(signed_bins)
  colors_tot = c(colors_neg, "#FFFFFF", "#FFFFFF", colors_pos)
  
  pheatmap(heat_data_w, scale="none", cellheight=10, cellwidth=10, file=out_file, color=colors_tot, breaks=seq(from=min_val, to=max_val, by=bin_size), cluster_cols=F)
  heat_data_w
}

getComparisons = function(qData){
  print(levels(qData$GROUP_ORIGINAL))
  ## read comparison matrix
  ## NS5A against control for BMS inhib, GSK inhib, No inhib respectively
  comparison = as.matrix(rbind(c(1,0,0,0,0,-1,0,0,0,0),c(0,1,0,0,0,0,-1,0,0,0),c(0,0,1,0,0,0,0,-1,0,0),c(0,0,0,1,0,0,0,0,-1,0),c(0,0,0,0,1,0,0,0,0,1),c(-1,1,0,0,0,0,0,0,0,0),c(-1,0,1,0,0,0,0,0,0,0),c(-1,0,0,1,0,0,0,0,0,0),c(-1,0,0,0,1,0,0,0,0,0))) 
  rownames(comparison)=c("BRG1D1-CTRLD1","BRG1D1-CTRLD2","BRG1D1-CTRLD3","BRG1D1-CTRLD5","BRG1D1-CTRLD10","BRG1D1-BRG1DALL","BRG1D2-BRG1DALL","BRG1D3-BRG1DALL","BRG1D5-BRG1DALL","BRG1D10-BRG1DALL")
  colnames(comparison)=levels(qData$GROUP_ORIGINAL)
  comparison
}

flattenTechRepeats = function(data){
  tmp1 = aggregate(Intensity ~ NormalizationGroup+Proteins+Sequence+Charge+IsotopeLabelType+Condition+BioReplicate, data=data, FUN=median, na.rm=T)
  run_rep = aggregate(cbind(Run,Raw.file) ~ BioReplicate, data = unique(data[,c("BioReplicate",'Run',"Raw.file")]), FUN=paste, collapse=",")
  tmp2 = merge(tmp1, run_rep, by="BioReplicate")
  tmp2[,c("Raw.file","Proteins","Sequence","Charge","Intensity","IsotopeLabelType","Condition","BioReplicate","Run","NormalizationGroup")]
}

######################################
## CALL ##############################

setwd("~/Projects/HPCKrogan/Data/HeartPPG/")

PDF=T
SAMPLE=F
SAMPLE_SIZE=100
NORMALIZE=T
PREPROCESS=T
MSSTATS=T
NORMALIZATION_METHOD="scale" ##scale, quantile, cyclicloess, toref
NORMALIZE_BATCH="TOREF" ##CONDITION, ALL, TOREF
FILL_MISSING=F
SPECIES="MOUSE"
FDR = 0.05 
LFC = 2
FLATTEN_TECH_REPEATS=F
REMOVE_PROTEINGROUPS=F
BACKGROUND_REF = "P99024|Q7TMM9|P68372|P99024|P68369|P68372|P99024|P68372|P99024|Q7TMM9|Q7TMM9|P68369" ## TUBULINS & GAPDH MOUSE
BAIT_REF = "Q3TKT4" ## BRG1 MOUSE
AGGREGATE_FUN=median

keys_file='data/input/031114-SH-22-69/031114-sh-22-69-BRG1-keys.txt'
data_file='data/input/031114-SH-22-69/052814-sh-22-69-evidence.txt'

contrast_file = 'data/input/031114-SH-22-69/BRG1_ALL_TOD0.txt'
BASE_NAME=sprintf("BRG1_ALLTOD0_N-%s_M-%s",NORMALIZE_BATCH, FILL_MISSING)

contrast_file = 'data/input/031114-SH-22-69/BRG1_MOCK_PER_D_contrast.txt'
BASE_NAME=sprintf("BRG1_MOCK_PER_D_N-%s_M-%s",NORMALIZE_BATCH, FILL_MISSING)

contrast_file = 'data/input/031114-SH-22-69/BRG1_MOCK_contrasts.txt'
BASE_NAME=sprintf("BRG1_MOCK_N-%s_M-%s",NORMALIZE_BATCH, FILL_MISSING)

PRINT_DIR=sprintf("Summary/20140618/%s",BASE_NAME)
dir.create(file.path(PRINT_DIR), showWarnings = FALSE)
out_file = sprintf('%s/%s.txt',PRINT_DIR,BASE_NAME)

main = function(){
  
  data = read.delim(data_file, stringsAsFactors=F)
  
  if(PREPROCESS){
    keys = read.delim(keys_file, stringsAsFactors=F)
    keys$Run = gsub('-|_','.',keys$Run)
    #rng_data = as.vector(keys$Run)
    
    data_f = maxqFilterAndJoin(keys=keys,data=data)
    if(PDF){
      pdf(file=gsub('.txt','-dist.pdf',out_file), width=7, height=7)
    }
    
    if(FLATTEN_TECH_REPEATS){
      data_ff = flattenTechRepeats(data_f)
    }else{
      data_ff = data_f
    }
    plotIntensityData(data_f)
    if(PDF){
      dev.off()
    }
    if(NORMALIZE){
      ## HERE!
      if(NORMALIZE_BATCH=="CONDITION"){
        data_n = normalizePerCondition(data_ff, NORMALIZATION_METHOD=NORMALIZATION_METHOD, FILL_MISSING)  
      }else if(NORMALIZE_BATCH=="ALL"){
        data_n = normalizeSingle(data_ff, NORMALIZATION_METHOD=NORMALIZATION_METHOD, FILL_MISSING)
      }else if(NORMALIZE_BATCH=="TOREF"){
        data_n = normalizeToProtein(data_ff, keys, bckgr_ref=BACKGROUND_REF, bait_ref=BAIT_REF, PDF=T, PRINT_DIR=PRINT_DIR)
      }
      out_file = gsub('.txt',sprintf('-%s.txt',NORMALIZATION_METHOD),out_file)
      if(PDF){
        pdf(file=gsub('.txt','-dist.pdf',out_file), width=7, height=7)
      }
      plotIntensityData(data_n)
      if(PDF){
        dev.off()
      }
      data_f=data_n
    }
    if(SAMPLE){
      data_f_sample = sampleProteins(data_f, SAMPLE_SIZE)
      data_f = data_f_sample
    }
    dmss = dataToMSSFormat(data_f)
    print(sqldf("select Condition, count(distinct(ProteinName)) as 'Proteins' from dmss group by Condition"))
    write.table(dmss,file=out_file, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
  }else{
    dmss = data
  }
  
  ## ABSTRACT FURTHER
  if(MSSTATS){
    dmss_miss = padMissing(dmss)
    qData = dataProcess(dmss_miss, normalization=F)  
    comparisons = as.matrix(read.delim(contrast_file, stringsAsFactors=F))
    results = groupComparison(contrast.matrix=comparisons, data=qData, labeled=F)  
    mss_out = annotate_with_uniprot(data=results$ComparisonResult, species=SPECIES, key="Protein")
    write.table(mss_out, file=gsub('.txt','-res.txt',out_file), eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
  }else{
    mss_out = read.delim(gsub('.txt','-res.txt',out_file), stringsAsFactors=F)
  }

  heat_data_w = plotHeat(mss_out, gsub('.txt','-res.pdf',out_file))
 
  sign_proteins = mss_out[(!is.na(mss_out$log2FC) & mss_out$adj_pvalue <= FDR & abs(mss_out$log2FC) >= LFC) | mss_out$Protein==BAIT_REF, 'Protein']
  mss_F = mss_out[mss_out$Protein %in% sign_proteins, ]
  heat_data_w = plotHeat(mss_F, gsub('.txt','-sign-res.pdf',out_file))
}

main()
