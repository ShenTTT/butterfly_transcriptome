# Purpose
These scripts are used to trim adaptor sequences from raw RNA-seq or sRNA-seq data.

# Usage
1. Install trimmomatic in your system and change the path in the scripts accordingly.

2. Enter the directory with your raw sequencing files. For pair-end RNA-seq, two fastq files: sample_1.fq.gz and sample_2.fq.gz are required per sample. For single-end sRNA-seq, only one fastq file is needed. The universal adaptor sequences (TruSeq3.fa) are provided by the trimmomatic package.

3. Run the script by typing: bash script.sh


