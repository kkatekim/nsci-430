#!/bin/bash

echo "Starting!"

python3 fix_expression_file.py

echo "Variants found!"

python3 do_analysis.py -i "../data/proteinDomain/RV_brain_high.csv" -n "RV_brainHigh" 
python3 do_analysis.py -i "../data/proteinDomain/RV_brain_medium.csv" -n "RV_brainMedium" 
python3 do_analysis.py -i "../data/proteinDomain/RV_brain_low.csv" -n "RV_brainLow"

echo "Brain counts done!"

python3 do_analysis.py -i "../data/proteinDomain/RV_all_high.csv" -n "RV_allHigh"
python3 do_analysis.py -i "../data/proteinDomain/RV_all_medium.csv" -n "RV_allMedium"
python3 do_analysis.py -i "../data/proteinDomain/RV_all_low.csv" -n "RV_allLow"

echo "All counts done!"

Rscript lr_protein_domain.R "RV_brainHigh"
Rscript lr_protein_domain.R "RV_brainMedium"
Rscript lr_protein_domain.R "RV_brainLow"

echo "Brain LR done!"

Rscript lr_protein_domain.R "RV_allHigh"
Rscript lr_protein_domain.R "RV_allMedium"
Rscript lr_protein_domain.R "RV_allLow"

echo "All LR done!"
