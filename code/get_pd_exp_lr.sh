#!/bin/bash

echo "Starting!"

python3 fix_expression_file.py -i "../data/proteinDomain/variants_in_protein_domain.csv" -r "pd_"

echo "Variants found!"

python3 do_analysis.py -i "../data/proteinDomain/pd_brain_high.csv" -n "RV_pd_brainHigh" 
python3 do_analysis.py -i "../data/proteinDomain/pd_brain_medium.csv" -n "RV_pd_brainMedium" 
python3 do_analysis.py -i "../data/proteinDomain/pd_brain_low.csv" -n "RV_pd_brainLow"

echo "Brain counts done!"

python3 do_analysis.py -i "../data/proteinDomain/pd_all_high.csv" -n "RV_pd_allHigh"
python3 do_analysis.py -i "../data/proteinDomain/pd_all_medium.csv" -n "RV_pd_allMedium"
python3 do_analysis.py -i "../data/proteinDomain/pd_all_low.csv" -n "RV_pd_allLow"

echo "All counts done!"

Rscript lr_protein_domain.R "RV_pd_brainHigh"
Rscript lr_protein_domain.R "RV_pd_brainMedium"
Rscript lr_protein_domain.R "RV_pd_brainLow"

echo "Brain LR done!"

Rscript lr_protein_domain.R "RV_pd_allHigh"
Rscript lr_protein_domain.R "RV_pd_allMedium"
Rscript lr_protein_domain.R "RV_pd_allLow"

echo "All LR done!"
