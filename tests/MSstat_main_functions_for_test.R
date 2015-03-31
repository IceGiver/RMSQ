#########################
## MAIN FUNCTIONS #######

filterData = function(data, config){
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
        if(!is.null(config$filters$modification)){
                cat(sprintf("\tMODIFICATIONS\t%s\n",config$filters$modification))
                if(config$filters$modification == 'UB'){
                        data_f = data_f[Modifications %like% 'GlyGly']
                }else if(config$filters$modification == 'PH'){
                        data_f = data_f[Modifications %like% 'Phospho']
                }
        }
        return(data_f)
}

aggregateData = function(data_w, keys, config, castFun){
        cat(">>AGGREGATING TECHNICAL REPEATS\n")
        ##  we need to convert to long format for more easy aggregation
        data_l = meltMaxQToLong(data_w, na.rm = T)
        keysagg = flattenKeysTechRepeats(keys)
        keysmapping = merge(keys, keysagg[,c('BioReplicate','RawFile')], by='BioReplicate')
        setnames(keysmapping,c('RawFile.x','RawFile.y'),c('RawFile','RawFileCombined'))
        data_l_combined = merge(data_l, keysmapping, by=c('RawFile','IsotopeLabelType'))
        aggregate_fun = tryCatch(eval(parse(text=config$aggregation$aggregate_fun)), 
                                 error = function(e) cat("\tWARNING: argument for aggregate_fun not a valid R function - defaulting to:\t 'max'"), 
                                 finally=max) 
        
        data_l_combined_agg = data.table(aggregate(Intensity ~ RawFileCombined + Proteins + Sequence + Charge + IsotopeLabelType, FUN=sum, data=data_l_combined))
        setnames(data_l_combined_agg,'RawFileCombined','RawFile')
        data_w_agg = castMaxQToWide(data_l_combined_agg)
        return(list(data_w_agg = data_w_agg, keys_agg = keysagg))
}

## returns data tabel in wide format
normalizeData = function(data_w, config){
        cat(">> NORMALIZING\n")
        
        if(grepl('scale|quantile|cyclicloess',config$normalization$method)){
                cat(sprintf("\tNORMALIZATION\t%s\n",config$normalization$method))
                data_fn = normalizeSingle(data_w=data_w, NORMALIZATION_METHOD=config$normalization$method)  
        }else if(grepl('reference',config$normalization$method) && !is.null(config$normalization$reference) && config$normalization$reference %in% unique(data_w$Proteins)){
                
                ## CURRENTLY UING MSSTATS METHODS FOR REFERENCE NORMALIZATION
                data_fn = data_w
                
                #     cat(sprintf("\tNORMALIZATION\tTO REFERENCE\t%s\n",config$normalization$reference))
                #     ref_files = keys[keys$NormalizationGroup == 'REFERENCE', 'RawFile']
                #     data_l_ref = data_l[data_l$Raw.file %in% ref_files & data_l$Intensity > 0 & is.finite(data_l$Intensity), ]
                #     data_l_nonref = data_l[!(data_l$Raw.file %in% ref_files) & data_l$Intensity > 0 & is.finite(data_l$Intensity), ]
                #     data_l_ref_n = normalizeToReference(data_l_ref=data_l_ref, ref_protein = config$normalization$reference, output_file = config$files$output)
                #     data_fn= rbind(data_l_ref_n, data_l_nonref)
                #     data_fn = castMaxQToWide(data_fn)
                
        }else{
                data_fn = data_w
        }
        
        return(data_fn)
}

runMSstats = function(dmss, contrasts, config){
        
        if(grepl('before', config$msstats$profilePlots)){
                mssquant = dataProcess(dmss, normalization=F, fillIncompleteRows = T)  
                dataProcessPlots(data=mssquant, type="ProfilePlot", featureName="Peptide", address=gsub('.txt','-before',config$files$output))
                dataProcessPlots(data=mssquant, type="QCPlot", address=gsub('.txt','-before',config$files$output))
        }
        
        if(!is.null(config$msstats$normalization_reference) & config$msstats$normalization_method == 'globalStandards'){
                normalization_refs = unlist(lapply(strsplit(config$msstats$normalization_reference, split = ','), FUN=trim))
                mssquant = dataProcess(dmss, normalization=config$msstats$normalization_method, nameStandards=normalization_refs , fillIncompleteRows=T)
        }else{
                mssquant = dataProcess(dmss, normalization=config$msstats$normalization_method , fillIncompleteRows=T)
        } 
        
        if(grepl('after', config$msstats$profilePlots)){
                dataProcessPlots(data=mssquant, type="ProfilePlot", featureName="Peptide", address=gsub('.txt','-after',config$files$output))
                dataProcessPlots(data=mssquant, type="QCPlot", address=gsub('.txt','-after',config$files$output))
        }
        
        if(!all(levels(mssquant$GROUP_ORIGINAL) == colnames(contrasts))){
                cat(sprintf('\tERROR IN CONTRAST COMPARISON: GROUP LEVELS DIFFERENT FROM CONTRASTS FILE\n\tGROUP LEVELS\t%s\n\tCONTRASTS FILE\t%s\n', paste(levels(mssquant$GROUP_ORIGINAL),collapse=','),paste(colnames(contrasts),collapse=',')))
                quit()
        } 
        
        cat(sprintf('\tFITTING CONTRASTS:\t%s\n',paste(rownames(contrasts),collapse=',')))
        write.table(mssquant, file=gsub('-results.txt','-mss-normalized.txt',config$files$output), eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
        results = groupComparison(data = mssquant, contrast.matrix = contrasts, labeled = as.logical(config$msstats$labeled), scopeOfBioReplication = config$msstats$scopeOfBioReplication, scopeOfTechReplication = config$msstats$scopeOfTechReplication, interference = as.logical(config$msstats$interference), equalFeatureVar = as.logical(config$msstats$equalFeatureVar), missing.action = config$msstats$missing_action)$ComparisonResult
        write.table(results, file=config$files$output, eol="\n", sep="\t", quote=F, row.names=F, col.names=T)  
        cat(sprintf(">> WRITTEN\t%s\n",config$files$output))
        return(results)
}


convertDataLongToMss = function(data_w, keys, config){
        cat(">> CONVERTING DATA TO MSSTATS FORMAT\n")
        data_l = meltMaxQToLong(data_w, na.rm = F)
        data_lk = mergeMaxQDataWithKeys(data_l, keys, by=c('RawFile','IsotopeLabelType'))
        dmss = dataToMSSFormat(data_lk)
        ## sanity check for zero's
        if(nrow(dmss[!is.na(dmss$Intensity) & dmss$Intensity == 0,]) > 0){
                dmss[!is.na(dmss$Intensity) & dmss$Intensity == 0,]$Intensity = NA
        } 
        write.table(dmss, file=gsub('.txt','-mss.txt',config$files$data), eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
        return(dmss)  
}

writeExtras = function(results, config){
        cat(">> ANNOTATING\n")
        if(config$output_extras$biomart){
                #mart = useMart(biomart = 'unimart')
                mart = useMart(biomart = 'unimart', dataset = 'uniprot', verbose = T)
                mart_anns = AnnotationDbi::select(mart, keytype='accession', columns=c('accession','name','protein_name','gene_name','ensembl_id'), keys=as.character(unique(results$Protein)))
                if(nrow(mart_anns) > 0){
                        mart_anns = aggregate(. ~ accession, data=mart_anns, FUN=function(x)paste(unique(x),collapse=','))
                        results_ann = merge(results, mart_anns, by.x='Protein', by.y='accession', all.x=T)
                        cat(sprintf('\tCHANGED OUTPUT FILE TO\t%s\n',config$files$output))
                        write.table(results_ann, file=gsub('.txt','-ann.txt',config$files$output), eol="\n", sep="\t", quote=F, row.names=F, col.names=T)
                        cat(sprintf(">> WRITTEN\t%s\n",config$files$output)) 
                }else{
                        results_ann = results
                }
        }else{
                results_ann = results
                config$files$output = config$output_extras$msstats_output
        }
        
        lfc_lower = as.numeric(unlist(strsplit(config$output_extras$LFC,split=" "))[1])
        lfc_upper = as.numeric(unlist(strsplit(config$output_extras$LFC,split=" "))[2])
        ## select subset of labels for heatmap and volcan plots
        selected_labels = config$output_extras$comparisons
        if(is.null(selected_labels) || selected_labels=='all') selected_labels='*'
        
        ## select data points  by LFC & FDR criterium in single condition and adding corresponding data points from the other conditions
        sign_hits = significantHits(results_ann,labels=selected_labels,LFC=c(lfc_lower,lfc_upper),FDR=config$output_extras$FDR)
        sign_labels = unique(sign_hits$Label)
        cat(sprintf("\tSELECTED HITS FOR PLOTS WITH LFC BETWEEN %s AND %s AT %s FDR:\t%s\n",lfc_lower, lfc_upper, config$output_extras$FDR, nrow(sign_hits)/length(sign_labels))) 
        
        ## REPRESENTING RESULTS AS HEATMAP
        if(config$output_extras$heatmap){
                ## plot heat map for all contrasts
                heat_labels = prettyPrintHeatmapLabels(uniprot_acs=sign_hits$Protein,uniprot_ids=sign_hits$name, gene_names=sign_hits$gene_name)
                heat_data_w = plotHeat(mss_F = sign_hits, out_file =  gsub('.txt','-sign.pdf',config$files$output), names=heat_labels, cluster_cols=config$output_extras$heatmap_cluster_cols, display = config$output_extras$heatmap_display)  
        }
        
        if(config$output_extras$volcano){
                file_name = gsub('.txt','-volcano.pdf',config$files$output)
                volcanoPlot(results_ann[grep(selected_labels,results_ann$Label),], lfc_upper, lfc_lower, FDR=config$output_extras$FDR, file_name=file_name)  
        }
}

trim <- function (x) gsub("^\\s+|\\s+$", "", x)
