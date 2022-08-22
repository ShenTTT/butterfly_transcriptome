# Purpose
These scripts are used to find alternative splicing events in the transcriptome and identify differentially spliced genes (DSGs) between conditions.

# Identify alternative splicing events

1. Install rMATS and create a sample.txt file to group all the RNA-seq samples from which to detect alternative splicing events, according to the rMTAS manual:

https://github.com/Xinglab/rmats-turbo/blob/v4.1.2/README.md

2. Detect alternative splicing events without comparisons using the command in rmats_statoff.txt

# Differential splicing analysis

The analysis is to detect differential splicing events between two conditions, by specifying two groups of samples (--b1 and --b2), using the command in rmats_compare.txt

# Filtering of the splicing events

All detected alternative splicing events and differential splicing events are subjected to further filtering and processing in R studio. Use AS.R to process total alternative splicing events/genes without comparisons, and generate overall splicing profiles (heatmap & PCA). Use DSG.R for differential splicing events/genes comparing different conditions.



