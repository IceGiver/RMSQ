RMSQ
====

dependencies
------------

CRAN
- getopt
- yaml
- reshape2
- ggplot
- pheatmap
- RColorBrewer
- data.table

BIOCONDUCTOR
- biomaRt
- limma
- MSstats

executing the main pipeline
---------------------------
chmod +x R/MSStats_main.R
Call R/MSStats_main.R --config 	YAMLCONFIG or 
Call R/MSStats_main.R -c 		YAMLCONFIG

example: 
1. unzip tests/LF/TIP47_evidence.txt.zip into tests/LF/TIP47_evidence.txt 
2. R/MSStats_main.R --config tests/LF/TIP47.yml

CONFIGURATION
-------------
about the yaml format for config files : http://www.yaml.org 
Every block enables and configures a conceptual part of the pipeline.

files :

  keys : A key file formatted as in tests/LF/TIP47_keys.txt that describes the specifics of each MS run, identified by the Raw.file column in the MaxQuant evidence file. Normalization group can be set to BACKGROUND or REFERENCE when we normalze against bait peptide levels. 

  data : A MaxQuant evidence file.

  contrasts : A file formatted as tests/LF/contrasts.txt that describes the experimental design. The clumn headers need to correspond to entries in the 'Condition' column of the keys file. The rownames describe the comparisons to be made. The coefficients -1,1,0 indicate which conditions need to be contrasted or ignored respectively. Be aware that the columns (conditions)should be in the same order as the R factors (by default alphabetical).   
  
  output : The output file name that will be written by MSsstats. This file name and location will be used as a base name for all the plots. 

  sequence_type : modified/unmodified. If modified is used then the Modified.sequence, which includes PTMs from the MaxQuant data file will be used as a unique peptide identifier. If this is set to 'unmodified' then the 'Sequence', column will be used as peptide identifier. This column does not contain PTMs. 

data:
  
  enabled : 0/1. Reads all data files. Can be disabled if you want to use the pipeline to produce plots only

filters: 
  
  enabled : 0/1. Enables this block


  contaminants : 0/1. Filters contaminants as defined by MaxQuant (REV__* and CON__*)

  protein_groups : remove/explode/ignore. Whether composite protein groups need to be filtered out (remove) or ignored. Currently explode is not implemented.

normalization :
  
  enabled : 0/1. Enables this block
  method : quantile/none/scale/reference. quantile/none/scale map to the normalizeBetweenArrays function of Limma. Reference performs normalization to a reference peptide's levels. Fill in reference with a Protein name for this to work. Typically used for pull downs where the bait level shoudl be the reference. 

  fill_missing : 0/1. Impute missing values or not. Currently the minimal observed intensity per run is used as the fill in value. 
  
  aggregate_tr : 0/1. Aggregate technical replicates into one 'virtual' run per biological replicate or not. Technical replicates are aggregated using the maximum observed value of a peptide across the technial repeats. 
  
  reference : A reference protein for normalization. Should be a value in the Protein column of the MaxQuant File. 

msstats :
  
  enabled : 0/1. Enables this block
  
  msstats_input : A file formatted according the msstats guidelines. If set this overrides the data computed in memory by the previous steps.   

output_extras :

  enabled : 0/1. Enables this block
  
  msstats_output : A file formatted according as msstats output. If set this overrides the data computed in memory by the previous steps. Typically used to produce multiple plots without having to re-run MSstats every time. Set this to the [FILE]-ann.txt format if biomart is bypassed (next param).
  
  biomart : 0/1. Gets all gene names and descriptions from Biomart for output puprposes. Leave this on in the pipeline if plots are enabled otherwise plots will fail.
  
  comparisons : all / any grep expression that returns a subset of the comparisons described in the contrasts file. Used to only plot subsets of comparisons in the heatmap
  
  LFC : lower / upper bound for Log2Fold changes that need to presented in the heat map. For example '-2 2' will select proteins that are < -2 or > 2 in any comparison. The respective Fold Change for Protein in the other comparisons will be plotted as well if it becomes selected. For pull down data the lower cutoff can be set to -Inf to select only Proteins that were significantly present in the bait. 
  
  FDR : p-value for FDR cutoff. The selection criteria for display on the heat map is a combination of both LFC and FDR.   
  
  heatmap : 0/1. Plot heat map or not.  
  
  heatmap_cluster_cols : 0/1. Cluster the columns of the hat map or not. If set to 0 the columns of the heat map will be displayed in the order of the R factors (by default alphabetical) 
  
  volcano : 0/1. Plot volcano plot with boundaries corresponding to the LFC pair and FDR or not.  
