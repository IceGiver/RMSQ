#! /usr/bin/Rscript --vanilla

suppressMessages(library(data.table))
suppressMessages(library(seqinr))
suppressMessages(library(stringr))
suppressMessages(require(bit64))
suppressMessages(library(getopt))
suppressMessages(library(reshape2))
suppressMessages(library(biomaRt))
suppressMessages(library(limma))

###############################
## FILE AND LIB LOADING #######

#########################
## CONFIG LOADING #######

ALLOWED_COMMANDS = c('concat','convert-silac','keys','convert-sites','annotate','results-wide','mapback-sites')

spec = matrix(c(
  'verbose', 'v', 2, "integer", "",
  'help'   , 'h', 0, "logical", "available arguments (this screen)",
  'command'  , 'c', 1, "character", sprintf("command to run. Currently supported commands: %s",paste(ALLOWED_COMMANDS,collapse=',')),
  'files'  , 'f', 1, "character", "files to feed to command. accepts regexp but needs to be quoted",
  'output'  , 'o', 1, "character", "Output file",
  'proteome'  , 'p', 1, "character", "Reference Proteome FASTA file",
  'mapping'  , 'm', 1, "character", "mapping file produced by convert-sites",
  'biomart_db'  , 'b', 1, "character", "name of biomart database to use for mapping, default=hsapiens_gene_ensembl",
  'mod_type','t', 1, "character", "Modification type: ub"),
  byrow=T, ncol=5)

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
  
  for(f in 1:length(files)){
    file = files[f]
    tmp = fread(file, stringsAsFactors=F, colClasses = c(Intensity='character'))
    tmp$Intensity = as.numeric(tmp$Intensity)
    #tmp$Intensity L = as.numeric(tmp[,'Intensity L',with=F])
    
    unique_files_current = unique(tmp[['Raw file']])
    if(!is.null(intersect(unique_files_current,unique_files)) && length(intersect(unique_files_current,unique_files))>0) cat(sprintf('\tWARNING DUPLICATE RAW FILE ENTRIES IN FILE %s:\t%s\n',file, paste(intersect(unique_files_current, unique_files),collapse=',')))
    select_colnames = grep('Raw\ file|Intensity|Proteins|Modifications|Sequence|Modified\ sequence|Charge|Protein\ group\ IDs|id|Retention\ time|Reverse|Contaminant',colnames(tmp), ignore.case = F)
    tmp = tmp[,select_colnames,with=F]
    res = rbind(res, tmp)
    unique_files = c(unique_files, unique_files_current)  
  }
  
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


MQutil.ProteinToSiteConversion <- function (maxq_file, ref_proteome_file, output_file, mod_type='ub') {
  
  if(mod_type=='ub'){
    maxq_mod_residue='K\\(gl\\)'
    mod_residue = 'K'
  }else if(mod_type=='ph'){
    maxq_mod_residue='(S|T|Y)\\(ph\\)'  
    mod_residue = 'S|T|Y'
  }
  
  ## read in ref. proteome
  ref_proteome = read.fasta(file = ref_proteome_file, 
                            seqtype = "AA", as.string = T,
                            set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, strip.desc = FALSE)
  
  ######################
  ## make mod-site index
  
  p_seqs = c()
  p_names = c()
  p_annots = c()
  
  for(e in ref_proteome){
    p_seqs = c(p_seqs, e[1])
    p_names = c(p_names, attr(e,'name'))
    p_annots = c(p_annots, attr(e,'Annot'))
  }
  
  ref_table = data.table(names=p_names, annots=p_annots, seqs=p_seqs)
  ref_table[,uniprot_ac:=gsub('([a-z,0-9,A-Z]+\\|{1})([A-Z,0-9,\\_]+)(\\|[A-Z,a-z,0-9,_]+)','\\2',names)]
  
  indices = lapply(ref_table$seqs, function(x) as.vector(str_locate_all(x,pattern=mod_residue)[[1]][,1]))
  ptm_sites = lapply(ref_table$seqs, function(x) as.vector(str_match_all(x,pattern=mod_residue)))
  lengths = unlist(lapply(indices, FUN = length))
  keys = rep(ref_table$uniprot_ac, lengths)
  protein_indices = data.table(uniprot_ac=keys, ptm_site=unlist(ptm_sites), res_index = unlist(indices))
  
  #################################
  ## map mod sites in data to index 
  
  ## read in maxq. data
  maxq_data = fread(maxq_file)
  unique_peptides_in_data = unique(maxq_data[,c('Proteins','Modified sequence'),with=F])
  setnames(unique_peptides_in_data,'Modified sequence','sequence')
  
  mod_sites = c()
  mod_seqs = c()
  
  for(i in 1:nrow(unique_peptides_in_data)){
    entry = unique_peptides_in_data[i,]
    peptide_seq = entry$sequence
    ## cleanup the sequence (removing all modifications) for matching the protein sequence
    peptide_seq_clean = gsub('[a-z,0-9,\\(,\\),_]','', peptide_seq)
    mod_sites_in_peptide = str_locate_all(string = peptide_seq, pattern = maxq_mod_residue)[[1]][,1]
    
    if(length(mod_sites_in_peptide)>0){
      uniprot_acs = entry$Proteins
      uniprot_acs = str_split(string = uniprot_acs, pattern = ';')[[1]]
      
      for(uac in uniprot_acs){
        protein_seq = ref_table[uniprot_ac==uac,]$seqs
        if(length(protein_seq)>0){
          ## get the position of the peptide in the protein sequence
          peptide_index_in_protein = str_locate(protein_seq, peptide_seq_clean)[[1]][1]
          
          for(m in 1:length(mod_sites_in_peptide)){
            mod_site = mod_sites_in_peptide[m]   
            peptide_seq_before_site = str_sub(peptide_seq, 1, mod_site-1)
            ## count all AA (not counting all modifications) before the modification to get the relative position of the modification in the peptide sequence
            residues_before_site = str_count(string = peptide_seq_before_site, pattern = 'A|C|D|E|F|G|H|I|K|L|M|N|P|Q|R|S|T|V|W|Y')
            mod_site_index_in_protein = peptide_index_in_protein+residues_before_site
            protein_mod_sites = protein_indices[uniprot_ac==uac,]
            if(mod_site_index_in_protein %in% protein_mod_sites$res_index){
              #cat(sprintf('%s\n',mod_site_id))
              mod_res = protein_mod_sites[res_index==mod_site_index_in_protein,ptm_site]
              mod_site_id = sprintf('%s_%s%s', uac, str_sub(protein_seq,mod_site_index_in_protein,mod_site_index_in_protein), mod_site_index_in_protein)
              mod_sites = c(mod_sites, mod_site_id)
              mod_seqs = c(mod_seqs, peptide_seq)
              stopifnot(length(mod_sites)==length(mod_seqs))
            }else{
              cat(sprintf('MISMATCH\t%s\n\tPEPTIDE_SEQ\t%s\n\tMOD_SITE\t%s\n\tPEPTIDE_IDX_IN_PROTEIN\t%s\n\tRESIDUES_BEFORE_SITE\t%s\n\tPROTEIN_SEQ\t%s\n',mod_site_id, peptide_seq, mod_site, peptide_index_in_protein, residues_before_site, protein_seq))
            } 
          }
        }
      }  
    }
  }
  
  mod_site_mapping = data.table(mod_sites, mod_seqs)
  mod_site_mapping_agg = aggregate(mod_sites ~ mod_seqs, mod_site_mapping, FUN=function(x)paste(x,collapse=','))
  
  setnames(maxq_data,'Modified sequence','mod_seqs')
  unmapped_mod_seqs = maxq_data[!(mod_seqs %in% mod_site_mapping_agg$mod_seqs) & grepl('(gl)',mod_seqs) & !grepl('REV__|CON__',Proteins),]
  unmapped_mod_seqs = unique(unmapped_mod_seqs[,c('mod_seqs','Proteins'),with=F])
  cat('UNABLE TO MAP\n')
  print(unmapped_mod_seqs)
  
  final_data = merge(maxq_data, mod_site_mapping_agg, by='mod_seqs')
  setnames(final_data,c('Proteins','mod_sites','mod_seqs'),c('Proteins_ref','Proteins','Modified sequence'))
  write.table(final_data, file = output_file, eol='\n', sep='\t',quote=F, row.names=F, col.names=T)
  
  ## write a mapping table
  protein_seq_mapping = unique(maxq_data[,c('Proteins','mod_seqs'),with=F])
  setnames(protein_seq_mapping,'Proteins','Protein')
  mapping_table = merge(protein_seq_mapping, mod_site_mapping_agg, by='mod_seqs', all=T)
  write.table(mapping_table, file=gsub('.txt','-mapping.txt',output_file), eol='\n', sep='\t',quote=F, row.names=F, col.names=T)
}


MQutil.annotate = function(input_file=opt$input, output_file=opt$output, uniprot_ac_col='Protein', group_sep=';', db='hsapiens_gene_ensembl'){
  
  cat(">> ANNOTATING\n")
  results = fread(input_file)
  mart = useMart(biomart = 'ensembl', dataset=db,verbose = T)
  ids = unique(results[,uniprot_ac_col,with=F])
  id_table = data.table(group_id=1:nrow(ids), uniprot_acs=ids)
  setnames(id_table, 2,uniprot_ac_col)
  results = merge(results, id_table, by=uniprot_ac_col)
  
  ids_split = c()
  uniprot_acs_split = c()
  
  for(id in 1:nrow(id_table)){
    id_row = id_table[id,]
    uniprot_ac_vec = unlist(str_split(unlist(id_row[,2,with=F]), pattern=group_sep))
    ids_split = c(ids_split,rep(unlist(id_row[,1,with=F]),times=length(uniprot_ac_vec)))
    uniprot_acs_split = c(uniprot_acs_split,uniprot_ac_vec)
  }
  
  id_table_split = unique(data.table(group_id=ids_split, uniprot_ac=uniprot_acs_split))
  
  mart_anns = getBM(mart = mart, attributes =c('uniprot_swissprot','uniprot_genename','description'), values=as.character(unique(id_table_split$uniprot_ac)), filter='uniprot_swissprot')
  mart_anns = aggregate(. ~ uniprot_swissprot, data=mart_anns, FUN=function(x)paste(unique(x),collapse=','))
  setnames(mart_anns, 'uniprot_swissprot', 'uniprot_ac')
  
  id_table_annotated = merge(id_table_split, mart_anns, 'uniprot_ac', all.x=T)
  id_table_annotated_flat = aggregate(. ~ group_id, data=id_table_annotated, FUN=function(x)paste(unique(x),collapse=','))
  results_out = merge(results, id_table_annotated_flat, by='group_id', all.x=T)
  
  unmapped = unique(results_out[is.na(results_out$uniprot_ac),uniprot_ac_col,with=F]) 
  cat(sprintf('UNMAPPED PROTEINS\t%s\n%s\n',nrow(unmapped),paste(unmapped$Protein,collapse=',')))
  
  write.table(results_out, file=output_file, sep='\t', quote=F, row.names=F, col.names=T)
  #return(results_out)
  
}

MQutil.resultsWide = function(input_file, output_file){
  input = fread(input_file)
  input_l = melt(data = input[,c('Protein', 'Label','log2FC','adj.pvalue'),with=F],id.vars=c('Protein', 'Label'))
  
  ## then cast to get combinations of LFCV/PVAl and Label as columns
  input_w = dcast.data.table( Protein ~ Label+variable, data=input_l, value.var=c('value'))
  write.table(input_w, file=output_file, eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
}

MQutil.mapSitesBack = function(input_file, mapping_file, output_file){
  input = fread(input_file)
  setnames(input,'Protein','mod_sites')
  mapping = fread(mapping_file)
  mapping = unique(mapping[!is.na(mod_sites),c('Proteins','mod_sites'),with=F])
  mapping = aggregate(Proteins ~ mod_sites, data=mapping, FUN=function(x)paste(x,collapse='T'))
  out = merge(input, mapping, by='mod_sites', all.x=T)
  write.table(out[,c(ncol(out),1:(ncol(out)-1)),with=F], file=output_file, eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
}

main <- function(opt){
  if(opt$command %in% ALLOWED_COMMANDS){
    cat(sprintf('>> EXECUTING:\t%s\n',opt$command))
    #loadLibs()
    if(opt$command == 'concat'){
      MQutil.concat(filenames=opt$files, output = opt$output)
    }else if(opt$command == 'convert-silac'){
      MQutil.SILACToLong(filename = opt$files, output = opt$output)
    }else if(opt$command == 'keys'){
      MQutil.getKeys(filename = opt$files, output = opt$output)
    }else if(opt$command == 'convert-sites'){
      MQutil.ProteinToSiteConversion (maxq_file = opt$files, output_file = opt$output, ref_proteome_file = opt$proteome, mod_type = opt$mod_type)
    }else if(opt$command == 'annotate'){
      MQutil.annotate(input_file = opt$files, output_file = opt$output , db = opt$biomart_db)
    }else if(opt$command == 'results-wide'){
      MQutil.resultsWide(input_file = opt$files, output_file = opt$output )
    }else if(opt$command == 'mapback-sites'){
      MQutil.mapSitesBack(input_file = opt$files, output_file = opt$output , mapping_file = opt$mapping)
    }
  }else{
    cat(sprintf('COMMAND NOT ALLOWED:\t%s\n',opt$command)) 
    cat(sprintf('ALLOWED COMMANDS:\t%s\n',paste(ALLOWED_COMMANDS,collapse=','))) 
  }
}

# opt$command = 'results-wide'
# opt$input = '~/Projects/HPCKrogan/Data/HIV-proteomics/results/20141124-ub-proteins/HIV-UB-SILAC-KROGAN-results.txt'
# opt$output = '~/Projects/HPCKrogan/Data/HIV-proteomics/results/20141124-ub-proteins/HIV-UB-SILAC-KROGAN-results.txt'

# opt$command = 'annotate'
# opt$input = '~/Projects/HPCKrogan/Data/HIV-proteomics/results/20141124-ub-proteins/HIV-UB-SILAC-KROGAN-results.txt'
# opt$output = '~/Projects/HPCKrogan/Data/HIV-proteomics/results/20141124-ub-proteins/HIV-UB-SILAC-KROGAN-results.txt'

# opt$command = 'convert-silac'
# opt$files = '~/Projects/HPCKrogan/Data/HIV-proteomics/Meena/abundance/HIV_vs_MOCK_PROTEIN_evidence.txt'
# opt$output = '~/Projects/HPCKrogan/Data/HIV-proteomics/Meena/abundance/HIV_vs_MOCK_PROTEIN_evidence_split.txt'

# opt$command = 'concat'
# opt$files = '~/Projects/HPCKrogan/Data/Mtb/Files/073113*evidence.txt'
# opt$output = '~/Projects/HPCKrogan/Data/HIV-proteomics/Meena/abundance/HIV_vs_MOCK_PROTEIN_evidence_split.txt'

# opt$command = 'convert-sites'
# opt$files =  '~/Projects/HPCKrogan/Data/HIV-proteomics/data//UB-silac//HIV-UB-SILAC-KROGAN-data.txt' 
# opt$output = '~/Projects/HPCKrogan/Data/HIV-proteomics/data//UB-silac//HIV-UB-SILAC-KROGAN-data-modK.txt'
# opt$mod_type = 'ub'
# opt$proteome = '~/Projects/HPCKrogan/Data/Uniprot/homo-sapiens-swissprot.fasta'

# opt$command = 'mapback-sites'
# opt$files =  '~/Projects/HPCKrogan/Data/HIV-proteomics/results//20141124-ub-sites/HIV-UB-SILAC-KROGAN-results.txt'
# opt$mapping =  '~/Projects/HPCKrogan/Data/HIV-proteomics/data/UB-silac/HIV-UB-SILAC-KROGAN-data-modK-mapping.txt'
# opt$output = '~/Projects/HPCKrogan/Data/HIV-proteomics/results//20141124-ub-sites/HIV-UB-SILAC-KROGAN-results.txt'

main(opt)
