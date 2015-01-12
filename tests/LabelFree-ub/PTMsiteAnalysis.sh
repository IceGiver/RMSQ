#!/usr/bin/sh
../../R/MaxQ_utilities.R -c convert-sites -f data/LabelFree-ub-evidence-sample.txt -o data/LabelFree-ub-evidence-sample-modK.txt -t ub -p ../files/homo-sapiens-swissprot.fasta
../../R/MSstats_main.R -c LabelFree-ub-test-sites.yaml 
../../R/MaxQ_utilities.R -c results-wide -f results/keepgroups-imputed/LabelFree-ub-evidence-results-modK.txt -o results/keepgroups-imputed/LabelFree-ub-evidence-results-modK-wide.txt
../../R/MaxQ_utilities.R -c mapback-sites -f results/keepgroups-imputed/LabelFree-ub-evidence-results-modK-wide.txt -o results/keepgroups-imputed/LabelFree-ub-evidence-results-modK-wide-mappedback.txt -m data/LabelFree-ub-evidence-sample-modK-mapping.txt
../../R/MaxQ_utilities.R -c annotate -f results/keepgroups-imputed/LabelFree-ub-evidence-results-modK-wide-mappedback.txt -o results/keepgroups-imputed/LabelFree-ub-evidence-results-modK-wide-mappedback-annotated.txt