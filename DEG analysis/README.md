# Purpose
These scripts are used to identify differentially expressed genes (DEGs) between conditions and generate overall gene expression profiles (heatmap & PCA).

# Gene quantification using Salmon

1. Install Salmon and generate Salmon indexed gentrome according to:
https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

2. Quantify transcript abundance using adaptor trimmed clean reads with salmon.sh. Please change the path before running the script.

# Differential gene expression analysis

The analysis is performed in R studio using the pipeline provided in DEG.R.


