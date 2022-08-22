#!/bin/bash 

for i in `ls -1 *_no_rRNA_no_tRNA_no_snRNA_no_snoRNA.fq | sed 's/_no_rRNA_no_tRNA_no_snRNA_no_snoRNA.fq//'` 

do 

echo "Select read length of 17-25 for ${i}"
awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 17 && length(seq) <= 25) {print header, seq, qheader, qseq}}' < ${i}_no_rRNA_no_tRNA_no_snRNA_no_snoRNA.fq > ${i}_CLEAN_17-25.fq


done

