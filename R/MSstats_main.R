#! /usr/bin/Rscript --vanilla
suppressMessages(library(getopt))
suppressMessages(library(yaml))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(org.Mm.eg.db))

#########################
## MIST MAIN FILE #######

spec = matrix(c(
  'verbose', 'v', 2, "integer", "",
  'help'   , 'h', 0, "logical", "available arguments (this screen)",
  'config'  , 'c', 1, "character", "configuration file in YAML format"),
  byrow=TRUE, ncol=5)

opt = getopt(spec = spec, opt = commandArgs(TRUE), command = get_Rscript_filename(), usage = FALSE, debug = FALSE)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}

# set source directory
args <- commandArgs(trailingOnly = F) 
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))

## load all externeal files
if(length(scriptPath)==0){
  source("R/MSstats_functions.R")
}else{
  source(paste(scriptPath,"/MSstats_functions.R",sep=""))  
}

main <- function(opt){
  config = tryCatch(yaml.load_file(opt$config), error = function(e) {cat(opt$config);break} )
  
  data = data.table(read.delim(config$files$data, stringsAsFactors=F))
  keys = read.delim(config$files$keys, stringsAsFactors=F)
  contrasts = as.matrix(read.delim(config$files$contrasts, stringsAsFactors=F))
  
  ## FILTERING
  if(config$filters$enabled){
    cat(">> FILTERING\n")
    if(config$filters$protein_groups == 'remove'){
      cat("\tPROTEIN GROUPS\tREMOVE\n")
      data_f = removeMaxQProteinGroups(data)  
    }else if(config$filters$protein_groups == 'explode'){
      cat("\tPROTEIN GROUPS\tEXPLODE\n")
      data_f = explodeMaxQProteinGroups(data)  
    }else{
      cat("\tPROTEIN GROUPS\tIGNORE\n")
      data_f = data
    }
    if(config$filters$contaminants){
      cat("\tCONTAMINANTS\tREMOVE\n")
      data_f = filterMaxqData(data_f)  
    } 
  }else{
    data_f = data
  }
  
  ## NORMALIZATION
  if(config$normalization$enabled){
    cat(">> NORMALIZING\n")
    data_w = castMaxQToWidePTM(data_f)
    if(config$normalization$fill_missing){
      cat("\tMISSING VALUES\tFILL\n")
      data_f = fillMissingMaxQEntries(data_w)
    } 
    
    if(grep('scale|quantile|cyclicloess',config$normalization$method)){
      cat(sprintf("\tNORMALIZATION\t%s\n",config$normalization$method))
      data_fn = normalizeSingle(data_w=data_f, NORMALIZATION_METHOD=config$normalization$method)  
    }else if(is.uniprot(config$normalization$method)){
      cat(sprintf("\tNORMALIZATION\tTO REFERENCE\t%s\n",config$normalization$method))
      cat('\tNOT YET IMPLEMENTED\t NO NORMALIZATION !!')
    }else{
      data_fn = data_f
    }
  }
  
  ## MSSTATS
  if(config$msstats$enabled){
    cat(">> MSSTATS\n")
    if(is.null(config$msstats$msstats_input)){
      data_l = meltMaxQToLong(data_fn)
      data_all = mergeMaxQDataWithKeys(data_l, keys, dataCol='Raw.file')
      dmss = dataToMSSFormat(data_all)  
    }else{
      cat(sprintf("\tREADING PREPROCESSED\t%s\n",config$msstats$msstats_input)) 
      dmss = read.delim(config$msstats$msstats_input, stringsAsfactors=F)
    }
    qData = dataProcess(dmss, normalization=F)    
    results = groupComparison(contrast.matrix=contrasts, data=qData, labeled=F)  
    write.table(results$ComparisonResult, file=config$files$output, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)  
    cat(sprintf(">> WRITTEN\t%s\n",config$files$output))
  }
  
  if(config$annotation$enabled){
    cat(">> ANNOTATING\n")
    db = switch(config$annotation$species, MOUSE=org.Mm.eg.db, HUMAN=org.Hs.eg.db)
    annotations = select(db, keys = as.character(unique(results$ComparisonResult$Protein)),columns=c("ENTREZID","SYMBOL","GENENAME"), keytype="UNIPROT")
    mss_out = merge(results$ComparisonResult,annotations, by.x='Protein', by.y='UNIPROT', all.x=T)
  }else{
    mss_out = results$ComparisonResult
  }
  write.table(mss_out, file=config$files$output, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)  
  mss_out = read.delim(file=config$files$output, stringsAsFactors=F)  
  cat(sprintf(">> WRITTEN\t%s\n",config$files$output))
  
  if(config$heatmap$enabled){
    cat(">> HEATMAP\n")
    if(!is.null(config$heatmap$msstats_output)){
      mss_out = read.delim(config$heatmap$msstats_output, stringsAsFactors=F)  
    }
    
    lfc_lower = as.numeric(unlist(strsplit(config$heatmap$LFC,split=" "))[1])
    lfc_upper = as.numeric(unlist(strsplit(config$heatmap$LFC,split=" "))[2])
    all_contrasts = significantHits(mss_out,labels='*',LFC=c(lfc_lower,lfc_upper),FDR=config$heatmap$FDR)
    cat(sprintf("\t SELECTED HITS BETWEEN %s AND %s AT %s FDR\t%s\n",lfc_lower, lfc_upper, config$heatmap$FDR, nrow(all_contrasts))) 
    heat_data_w = plotHeat(all_contrasts, gsub('.txt','-sign.pdf',config$files$output), names=paste(all_contrasts$Protein,all_contrasts$SYMBOL, sep=' | '), cluster_cols = config$heatmap$cluster_cols)
  } 
}

if(!exists('DEBUG') || DEBUG==F){
  opt = c(opt, config='tests/MSstats_main_test.yml')
  main(opt)
}else{
  main(opt)
} 

# #########
# ## CONFIG 
# 
# SAMPLE=F
# SAMPLE_SIZE=100
# NORMALIZE=T
# PREPROCESS=T
# MSSTATS=T
# NORMALIZATION_METHOD="quantile" ##scale, quantile, cyclicloess, toref
# NORMALIZE_BATCH="ALL" ##CONDITION, ALL, TOREF
# FILL_MISSING=T
# SPECIES="MOUSE"
# FDR = 0.05 
# LFC = 2
# REMOVE_PROTEINGROUPS=F
# VERSION=1
# 
# keys_file='data/input/031114-SH-22-69/031114-sh-22-69-BRG1-keys.txt'
# data_file='data/input/031114-SH-22-69/052814-sh-22-69-evidence.txt'
# contrast_file = 'data/input/031114-SH-22-69/BRG1_MOCK_PER_D_contrast.txt'
# 
# BASE_NAME=sprintf("BRG1_MOCK_N-%s_M-%s",NORMALIZE_BATCH, FILL_MISSING)
# PRINT_DIR=sprintf("summary/20140714/%s",BASE_NAME)
# dir.create(file.path(PRINT_DIR), showWarnings = FALSE)
# out_file = sprintf('%s/%s_v%s.txt',PRINT_DIR,BASE_NAME,VERSION)
# 
# #######
# ## MAIN
# 
# 
# 
# ###
# 
# mss_out = read.delim(gsub('.txt','-res.txt',out_file), stringsAsFactors=F)
# 
# 
# mss_out_with_genes = merge(mss_out, gene_mapping, by.x='Protein', by.y='UNIPROT', all.x=T)
# mss_out_with_genes = mss_out_with_genes[,c("Protein","Label","log2FC","adj_pvalue","uniprot_id","ENTREZID","SYMBOL","GENENAME","description")]
