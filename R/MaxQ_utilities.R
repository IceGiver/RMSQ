#! /usr/bin/Rscript --vanilla

require(bit64)

###############################
## FILE AND LIB LOADING #######

suppressMessages(library(getopt))

#########################
## CONFIG LOADING #######

ALLOWED_COMMANDS = c('concat','convert-silac','keys')

spec = matrix(c(
  'verbose', 'v', 2, "integer", "",
  'help'   , 'h', 0, "logical", "available arguments (this screen)",
  'command'  , 'c', 1, "character", sprintf("command to run. Currently supported commands: %s",paste(ALLOWED_COMMANDS,collapse=',')),
  'files'  , 'f', 1, "character", "files to feed to command. accepts regexp but needs to be quoted",
  'output'  , 'o', 1, "character", "Output file"),
  byrow=TRUE, ncol=5)

opt = getopt(spec = spec, opt = commandArgs(TRUE), command = get_Rscript_filename(), usage = FALSE, debug = FALSE)

# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE));
  q(status=1);
}

loadLibs = function(){
  cat(">> LOADING EXTERNAL FILES AND LIBRARIES\n")
  # set source directory
  args <- commandArgs(trailingOnly = F) 
  scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
  
  ## load all external files
  if(length(scriptPath)==0){
    source("R/MSstats_functions.R")
  }else{
    source(paste(scriptPath,"/MSstats_functions.R",sep=""))  
  }
}

MQutil.SILACToLong = function(filename, output){
  library(data.table)
  library(reshape2)
  file = Sys.glob(filename)
  cat(sprintf('\tPROCESSING:\n\t%s\n',paste(file,collapse='\n\t')))
  tmp = fread(file)
  tmp_long = reshape2::melt(tmp, measure.vars = c('Intensity L','Intensity H'))
  tmp_long[,Intensity:=NULL]
  setnames(tmp_long,'value','Intensity')
  setnames(tmp_long,'variable','IsotopeLabelType')
  setnames(tmp_long,'Raw file','Raw.file')
  levels(tmp_long$IsotopeLabelType) = c('L','H')
  tmp_long[!is.na(tmp_long$Intensity) && tmp_long$Intensity<1,]$Intensity=NA
  write.table(tmp_long, file=output, sep='\t', quote=F, row.names=F, col.names=T)
}

MQutil.concat = function(filenames, output){
  library(data.table)
  files = Sys.glob(filenames)
  cat(sprintf('\tPROCESSING:\n\t%s\n',paste(files,collapse='\n\t')))
  
  res = NULL
  unique_files = c()
  
  for(file in files){
    tmp = fread(file, stringsAsFactors=F, colClasses = c(Intensity='character'))
    tmp$Intensity = as.numeric(tmp$Intensity)
    #tmp$Intensity L = as.numeric(tmp[,'Intensity L',with=F])
    
    unique_files_current = unique(tmp[['Raw file']])
    if(!is.null(intersect(unique_files_current,unique_files)) && length(intersect(unique_files_current,unique_files))>0) cat(sprintf('\tWARNING DUPLICATE RAW FILE ENTRIES IN FILE %s:\t%s\n',file, paste(intersect(unique_files_current, unique_files),collapse=',')))
    res = rbind(res, tmp)
    unique_files = c(unique_files, unique_files_current)  
  }
  select_colnames = grep('Raw\ file|Intensity|Proteins|Modifications|Sequence|Modified\ sequence|Charge|Protein\ group\ IDs|id|Retention\ time|Reverse|Contaminant',colnames(res), ignore.case = F)
  res = res[,select_colnames,with=F]
  cat(sprintf('\tWRITING\t%s\n',output))
  write.table(res, file=output, eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
  cat(sprintf('\tWRITING\t%s\n',gsub('.txt','-keys.txt',output)))
  #print(unique_files)
  write.table(unique_files, file=gsub('.txt','-keys.txt',output), eol='\n', sep='\t', quote=F, row.names=F, col.names=F)
}

MQutil.getKeys = function(filename, output){
  library(data.table)
  file = Sys.glob(filename)
  cat(sprintf('\tPROCESSING:\n\t%s\n',paste(file,collapse='\n\t')))
  #tmp = data.table(read.delim(file, stringsAsFactors=F))
  tmp = fread(file, stringsAsFactors=F)
  write.table(unique(tmp$Raw.file),file=output, eol='\n', sep='\t', quote=F, row.names=F, col.names=F)
  cat(sprintf('\tWRITTEN\t%s\n',output))
}

main <- function(opt){
  if(opt$command %in% ALLOWED_COMMANDS){
    cat(sprintf('>> EXECUTING:\t%s\n',opt$command))
    #loadLibs()
    if(opt$command == 'concat'){
      MQutil.concat(opt$files, opt$output)
    }else if(opt$command == 'convert-silac'){
      MQutil.SILACToLong(filename = opt$files, output = opt$output)
    }else if(opt$command == 'keys'){
      MQutil.getKeys(filename = opt$files, output = opt$output)
    }  
  }else{
    cat(sprintf('COMMAND NOT ALLOWED:\t%s\n',opt$command)) 
    cat(sprintf('ALLOWED COMMANDS:\t%s\n',paste(ALLOWED_COMMANDS,collapse=','))) 
  }
}

# opt$command = 'convert-silac'
# opt$files = '~/Projects/HPCKrogan/Data/HIV-proteomics/Meena/abundance/HIV_vs_MOCK_PROTEIN_evidence.txt'
# opt$output = '~/Projects/HPCKrogan/Data/HIV-proteomics/Meena/abundance/HIV_vs_MOCK_PROTEIN_evidence_split.txt'

# opt$command = 'concat'
# opt$files = '~/Projects/HPCKrogan/Data/Mtb/Files/073113*evidence.txt'
# opt$output = '~/Projects/HPCKrogan/Data/HIV-proteomics/Meena/abundance/HIV_vs_MOCK_PROTEIN_evidence_split.txt'

main(opt)
