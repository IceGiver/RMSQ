#! /usr/bin/Rscript --vanilla --default-packages=utils
suppressMessages(library(grDevices))
suppressMessages(library(stats))
suppressMessages(library(graphics))
suppressMessages(library(reshape2))
suppressMessages(library(limma))
suppressMessages(library(ggplot2))
#library(MSstats.daily)
suppressMessages(library(MSstats))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(data.table))

theme_set(theme_bw(base_size = 15, base_family="Helvetica"))

############
## FUNCTIONS

is.uniprotAc = function(identifier){
  grepl('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}',identifier)
}

filterMaxqData = function(data){
  data_selected = data[grep("CON__|REV__",data$Proteins, invert=T),]
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

cleanMissingMaxQEntries = function(data_l){
  data_l[data_l$Intensity <= 0 | is.infinite(data_l$Intensity) | is.nan(data_l$Intensity),]$Intensity=NA
  return(data_l)
}

flattenMaxQTechRepeats = function(data_l){
  
}

flattenKeysTechRepeats = function(keys){
  aggregate(. ~ BioReplicate, data=keys, FUN=function(x)paste(unique(x),collapse='|'))
}

mergeMaxQDataWithKeys = function(data_l, keys, dataCol='Raw.file'){
  setnames(data_l, dataCol, 'RawFile')
  unique_data = unique(data_l$RawFile)
  unique_keys = unique(keys$RawFile)
  keys_not_found = setdiff(unique_keys, unique_data)
  data_not_found = setdiff(unique_data, unique_keys)
  cat(sprintf("keys found: %s \t keys not in data file:\n%s\n", length(unique_keys)-length(keys_not_found), paste(keys_not_found,collapse='\t')))
  cat(sprintf("data found: %s \t data not in keys file:\n%s\n", length(unique_data)-length(data_not_found), paste(data_not_found, collapse='\t')))
  
  ## select only required attributes from MQ format
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

plotHeat = function(mss_F, out_file, labelOrder=NULL, names='Protein', cluster_cols=F){
  heat_data = data.frame(mss_F, names=names)
  #heat_data = mss_F[,c('uniprot_id','Label','log2FC')]
  heat_data_w = dcast(names ~ Label, data=heat_data, value.var='log2FC')
  #gene_names = uniprot_to_gene_replace(uniprot_ac=heat_data_w$Protein)
  rownames(heat_data_w) = heat_data_w$names
  heat_data_w = heat_data_w[,-1]
  heat_data_w[is.na(heat_data_w)]=0
  max_val = ceiling(max(heat_data_w))
  min_val = floor(min(heat_data_w))
  extreme_val = max(max_val, abs(min_val))
  if(extreme_val %% 2 != 0) extreme_val=extreme_val+1
  bin_size=2
  signed_bins = (extreme_val/bin_size)
  colors_neg = rev(colorRampPalette(brewer.pal("Blues",n=extreme_val/bin_size))(signed_bins))
  colors_pos = colorRampPalette(brewer.pal("Reds",n=extreme_val/bin_size))(signed_bins)
  colors_tot = c(colors_neg, colors_pos)
  
  if(is.null(labelOrder)){
    pheatmap(heat_data_w, scale="none", cellheight=10, cellwidth=10, file=out_file, color=colors_tot, breaks=seq(from=-extreme_val, to=extreme_val, by=bin_size), cluster_cols=cluster_cols, fontfamily="mono")  
  }else{
    heat_data_w = heat_data_w[,labelOrder]
    pheatmap(heat_data_w, scale="none", cellheight=10, cellwidth=10, file=out_file, color=colors_tot, breaks=seq(from=-extreme_val, to=extreme_val, by=bin_size), cluster_cols=cluster_cols, fontfamily="mono")
  }
  
  heat_data_w
}

significantHits = function(mss_results, labels='*', LFC=c(-2,2), FDR=0.05){
  ## get subset based on labels
  selected_results = mss_results[grep(labels,mss_results$Label), ]
  cat(sprintf('\tAVAILABLE LABELS FOR HEATMAP:\t%s\n',paste(unique(mss_results$Label), collapse=',')))
  cat(sprintf('\tSELECTED LABELS FOR HEATMAP:\t%s\n',paste(unique(selected_results$Label), collapse=',')))
  significant_proteins = selected_results[(!is.na(selected_results$log2FC) & selected_results$adj.pvalue <= FDR & (selected_results$log2FC >= LFC[2] | selected_results$log2FC <= LFC[1])) , 'Protein']
  significant_results = selected_results[selected_results$Protein %in% significant_proteins, ]
  return(significant_results)
}

logScale = function(data, format='wide', base=2){
  if(format=='wide'){
    data[,4:ncol(data)] = log(data[,4:ncol(data)], base=base)
  }else{
    data$Intensity = log(data$Intensity, base=base)
  }
  return(data)
}

###################################
## doesnt work so far with new code


peptideDistribution = function(data_l, output_file, PDF=T){
  ## look at peptide distribution of all proteins
  if(PDF) pdf(output_file, width=10, height=7)
  p = ggplot(data=data_l, aes(x=Raw.file, y=Intensity))
  print(p + geom_boxplot() +
          theme(axis.text.x=element_text(angle=-90)))  
  if(PDF) dev.off()
}

peptideIntensityPerFile = function(ref_peptides, output_file, PDF=T){
  if(PDF) pdf(output_file, width=10, height=7)
  p = ggplot(data=ref_peptides, aes(x=Raw.file, y=Intensity, group=protein_id))
  print(p + geom_line(colour="grey") + 
          stat_summary(aes(group=1), geom="point", fun.y=median, shape=17, size=3, colour="darkred") + 
          ylab("Intensity") +
          theme(axis.text.x=element_text(angle=-90)))
  if(PDF) dev.off()
}

volcanoPlot = function(mss_results_sel, lfc_upper, lfc_lower, FDR, file_name='', PDF=T){
  
  min_x = -ceiling(max(abs(mss_results_sel$log2FC)))
  max_x = ceiling(max(abs(mss_results_sel$log2FC)))
  min_pval = min(mss_results_sel[mss_results_sel$bb > 0,]$adj.pvalue)
  mss_results_sel[mss_results_sel$adj.pvalue == 0,]$adj.pvalue = min_pval
  max_y = ceiling(-log2(min(mss_results_sel[mss_results_sel$adj.pvalue > 0,]$adj.pvalue)))
  
  # top
  
  if(PDF) pdf(file_name, width=7, height=7)
  p = ggplot(mss_results_sel, aes(x=log2FC,y=-log2(adj.pvalue)))
  print(p + geom_point(colour='grey') + 
    geom_point(data = mss_results_sel[mss_results_sel$adj.pvalue <= FDR & mss_results_sel$log2FC>=lfc_upper,], aes(x=log2FC,y=-log2(adj.pvalue)), colour='red', size=2) +
    geom_point(data = mss_results_sel[mss_results_sel$adj.pvalue <= FDR & mss_results_sel$log2FC<=lfc_lower,], aes(x=log2FC,y=-log2(adj.pvalue)), colour='blue', size=2) +
    geom_vline(xintercept=c(lfc_lower,lfc_upper), lty='dashed') + 
    geom_hline(yintercept=-log2(FDR), lty='dashed') + 
    xlim(min_x,max_x) + 
    ylim(0,max_y))
  if(PDF) dev.off()
}

prettyPrintHeatmapLabels = function(uniprot_acs, uniprot_ids, gene_names){
  uniprot_ids_trunc = gsub('([A-Z,0-9]+)_([A-Z,0-9]+)','\\1',uniprot_ids)
  longest_id = max(nchar(uniprot_ids_trunc))
  tmp_frame = data.frame(t=uniprot_ids_trunc, s=longest_id-nchar(uniprot_ids_trunc)+1, g=gene_names,a=uniprot_acs, stringsAsFactors=F)
  tmp_frame[is.na(tmp_frame$t),]$t=tmp_frame[is.na(tmp_frame$t),]$a
  result = apply(tmp_frame, 1, function(x)paste0(x[1],paste(rep(' ',x[2]),collapse=''),x[3]))
  return(result)
}

normalizeToReference = function(data_l_ref, ref_protein, PDF=T, output_file){
  
  data_l_ref$Intensity = log2(data_l_ref$Intensity)  
  
  ## SAMPLES WITH REF PROTEINS BEFORE 
  peptideDistribution(data_l_ref, gsub('.txt','-protein-dist-b.pdf',output_file))
  
  ## compute the average value for background reference peptides over all samples (background+bait)
  ## normalize against background reference proteins
  ref_peptides = data_l_ref[grep(ref_protein, data_l_ref$Proteins),]
  ref_peptides = data.frame(ref_peptides, protein_id=as.character(paste0(ref_peptides$Sequence,ref_peptides$Charge)), stringsAsFactors=F)
  ref_peptides_counts = aggregate(Raw.file ~ protein_id,data=ref_peptides,FUN=function(x)length(unique(x)))
  unique_files = length(unique(data_l_ref$Raw.file))
  complete_observations = ref_peptides_counts[ref_peptides_counts$Raw.file == unique_files,'protein_id']
  ref_peptides_complete = ref_peptides[ref_peptides$protein_id %in% complete_observations,]
  
  ## look at their individual signal over replicates
  peptideIntensityPerFile(ref_peptides_complete[!is.na(ref_peptides_complete$Intensity) & is.finite(ref_peptides_complete$Intensity),], gsub('.txt','-peptide-signal-b.pdf',output_file))
  
  ref_peptides_avg_per_rep = aggregate(Intensity ~ Raw.file, data=ref_peptides_complete, FUN=median)
  ref_peptides_avg_all = median(ref_peptides_avg_per_rep$Intensity)
  ref_peptides_avg_per_rep = data.frame(Raw.file=ref_peptides_avg_per_rep$Raw.file, correction=ref_peptides_avg_all-ref_peptides_avg_per_rep$Intensity)
  ref_peptides_complete = merge(ref_peptides_complete, ref_peptides_avg_per_rep, by='Raw.file', all.x=T)
  ref_peptides_complete = data.frame(ref_peptides_complete, IntensityCorrected=ref_peptides_complete$Intensity+ref_peptides_complete$correction)
  
  ref_peptides_complete$Intensity=ref_peptides_complete$IntensityCorrected
  peptideIntensityPerFile(ref_peptides_complete[!is.na(ref_peptides_complete$Intensity) & is.finite(ref_peptides_complete$Intensity),], gsub('.txt','-peptide-signal-a.pdf',output_file))

  data_l_ref = merge(data_l_ref, ref_peptides_avg_per_rep, by='Raw.file', all.x=T)
  data_l_ref$Intensity=data_l_ref$Intensity+data_l_ref$correction
  data_l_ref = data_l_ref[,-(ncol(data_l_ref))]
  
  ## SAMPLES WITH REF PROTEINS AFTER
  peptideDistribution(data_l_ref, gsub('.txt','-protein-dist-a.pdf',output_file))
  data_l_ref$Intensity=2^(data_l_ref$Intensity)
  
  return(data_l_ref)
}
