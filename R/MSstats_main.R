#! /usr/bin/Rscript --vanilla

###############################
## FILE AND LIB LOADING #######

cat(">> LOADING EXTERNAL FILES AND LIBRARIES\n")

suppressMessages(library(getopt))
suppressMessages(library(yaml))
suppressMessages(library(biomaRt))

# set source directory
args <- commandArgs(trailingOnly = F) 
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))

## load all external files
if(length(scriptPath)==0){
  source("R/MSstats_functions.R")
}else{
  source(paste(scriptPath,"/MSstats_functions.R",sep=""))  
}

#########################
## CONFIG LOADING #######

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

###############
## MAIN #######

main <- function(opt){
  cat(">> MSSTATS PIPELINE\n")
  config = tryCatch(yaml.load_file(opt$config), error = function(e) {cat(opt$config);break} )
  
  if(config$data$enabled){
    cat(">> LOADING DATA\n")
    data = data.table(read.delim(config$files$data, stringsAsFactors=F))
    keys = read.delim(config$files$keys, stringsAsFactors=F)
    contrasts = as.matrix(read.delim(config$files$contrasts, stringsAsFactors=F))  
  }
  
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
    
    ## formatting
    if(config$files$sequence_type == 'modified'){
      data_w = castMaxQToWidePTM(data_f)  
    }else if(config$files$sequence_type == 'unmodified'){
      data_w = castMaxQToWide(data_f)  
    }
    data_l = meltMaxQToLong(data_w)
    
    ## flattening tech repeats per biological replicate
    if(config$normalization$aggregate_tr){
      cat("\tAGGREGATING TECHNICAL REPEATS\n")
      newkeys = flattenKeysTechRepeats(keys)
      keys_tmp = merge(keys, newkeys[,c('BioReplicate','RawFile')], by='BioReplicate')
      setnames(keys_tmp,c('RawFile.x','RawFile.y'),c('RawFile','RawFileCombined'))
      tmp = merge(data_l, keys_tmp, by.x='Raw.file',by.y='RawFile')
      tmp = aggregate(Intensity ~ RawFileCombined + Proteins + Sequence + Charge, FUN=max, data=tmp)
      data_l = tmp
      setnames(data_l,'RawFileCombined','Raw.file')
      data_w = castMaxQToWide(data_l)
      keys = newkeys
    }
    
    ## normalization
    if(grepl('scale|quantile|cyclicloess',config$normalization$method)){
      cat(sprintf("\tNORMALIZATION\t%s\n",config$normalization$method))
      data_fn = normalizeSingle(data_w=data_w, NORMALIZATION_METHOD=config$normalization$method)  
    }else if(grepl('reference',config$normalization$method) & config$normalization$reference %in% data_l$Proteins){
      cat(sprintf("\tNORMALIZATION\tTO REFERENCE\t%s\n",config$normalization$reference))
      ref_files = keys[keys$NormalizationGroup == 'REFERENCE', 'RawFile']
      data_l_ref = data_l[data_l$Raw.file %in% ref_files & data_l$Intensity > 0 & is.finite(data_l$Intensity), ]
      data_l_nonref = data_l[!(data_l$Raw.file %in% ref_files) & data_l$Intensity > 0 & is.finite(data_l$Intensity), ]
      data_l_ref_n = normalizeToReference(data_l_ref=data_l_ref, ref_protein = config$normalization$reference, output_file = config$files$output)
      data_fn= rbind(data_l_ref_n, data_l_nonref)
      data_fn = castMaxQToWide(data_fn)
    }else{
      data_fn = data_w
    }
    ## herafter we assume data_fn is in wide format
    
    if(config$normalization$fill_missing){
      cat("\tMISSING VALUES\tFILL\n")
      data_fn = fillMissingMaxQEntries(data_w)
    }
  }
  
  ## MSSTATS
  if(config$msstats$enabled){
    cat(">> MSSTATS\n")
    if(is.null(config$msstats$msstats_input)){
      data_l = meltMaxQToLong(data_fn)
      
      ## if missing entries were not filled we need to remove the 0, Inf and NaN values to replace them with NA
      if(!config$normalization$fill_missing) data_l = cleanMissingMaxQEntries(data_l)
      
      data_all = mergeMaxQDataWithKeys(data_l, keys, dataCol='Raw.file')
      dmss = dataToMSSFormat(data_all)  
    }else{
      cat(sprintf("\tREADING PREPROCESSED\t%s\n",config$msstats$msstats_input)) 
      dmss = read.delim(config$msstats$msstats_input, stringsAsfactors=F)
    }
    qData = dataProcess(dmss, normalization=F)    
    results = groupComparison(contrast.matrix=contrasts, data=qData, labeled=F)$ComparisonResult  
    write.table(results, file=config$files$output, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)  
    cat(sprintf(">> WRITTEN\t%s\n",config$files$output))
    mss_out = results
  }
  
  ## ANNOTATING RESULT FILE
  if(config$output_extras$enabled){
    cat(">> ANNOTATING\n")
    if(! is.null(config$output_extras$msstats_output)) results = read.delim(config$output_extras$msstats_output, stringsAsFactors=F)
    mart = useMart(biomart = 'unimart', dataset = 'uniprot')
    mart_anns = AnnotationDbi::select(mart, keytype='accession', columns=c('accession','name','protein_name','gene_name','ensembl_id'), keys=as.character(unique(results$Protein)))
    mart_anns = aggregate(. ~ accession, data=mart_anns, FUN=function(x)paste(unique(x),collapse=',')) 
    mss_out = merge(results, mart_anns, by.x='Protein', by.y='accession', all.x=T) 
    config$files$output = gsub('.txt','-ann.txt',config$files$output)
    cat(sprintf('\tCHANGED OUTPUT FILE TO\t%s\n',config$files$output))
    write.table(mss_out, file=config$files$output, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)  
    cat(sprintf(">> WRITTEN\t%s\n",config$files$output))
    
    lfc_lower = as.numeric(unlist(strsplit(config$output_extras$LFC,split=" "))[1])
    lfc_upper = as.numeric(unlist(strsplit(config$output_extras$LFC,split=" "))[2])
    
    ## REPRESENTING RESULTS AS HEATMAP
    if(config$output_extras$heatmap){
      
      ## select subset of labels for heatmap and volcan plots
      selected_labels = config$output_extras$comparisons
      if(is.null(selected_labels) || selected_labels=='all') selected_labels='*'
      
      ## select data points for heatmap by LFC & FDR criterium in single condition and adding corresponding data points from the other conditions
      sign_hits = significantHits(mss_out,labels=selected_labels,LFC=c(lfc_lower,lfc_upper),FDR=config$output_extras$FDR)
      sign_labels = unique(sign_hits$Label)
      cat(sprintf("\tSELECTED HITS FOR PLOTS WITH LFC BETWEEN %s AND %s AT %s FDR:\t%s\n",lfc_lower, lfc_upper, config$output_extras$FDR, nrow(sign_hits)/length(sign_labels))) 
      
      ## plot heat map for all contrasts
      heat_labels = prettyPrintHeatmapLabels(uniprot_acs=sign_hits$Protein,uniprot_ids=sign_hits$name, gene_names=sign_hits$gene_name)
      heat_data_w = plotHeat(sign_hits, gsub('.txt','-sign.pdf',config$files$output), names=heat_labels, cluster_cols=config$output_extras$heatmap_cluster_cols)  
    }
    
    if(config$output_extras$volcano){
      ## make a volcano plat per contrast
#       for(l in sign_labels){
#         mss_results_sel = mss_out[mss_out$Label == l, ]
#         file_name = gsub('.txt',sprintf('-%s.pdf',l),config$files$output)
#         volcanoPlot(mss_results_sel, lfc_upper, lfc_lower, FDR=config$output_extras$FDR, file_name =file_name)
#       }  

      file_name = gsub('.txt','-volcano.pdf',config$files$output)
      volcanoPlot(mss_out, lfc_upper, lfc_lower, FDR=config$output_extras$FDR, file_name=file_name)  
    }
  }
  
}

if(!exists('DEBUG') || DEBUG==F){
  opt = c(opt, config='tests/MSstats_main_test.yml')
  main(opt)
}else{
  main(opt)
} 
