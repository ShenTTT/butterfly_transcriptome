#!/bin/bash 

for i in `ls -1 *_Aligned.out.sam | sed 's/_Aligned.out.sam//'` 

do 

perl unitas_1.7.7.pl -i ${i}_Aligned.out.sam -s x -refseq miRNA_hairpin_prefix_wing.fa -refseq rRNA_prefix.fa -refseq tRNA_prefix.fa -refseq snRNA_prefix.fa -refseq snoRNA_prefix.fa -refseq RNA_prefix.fa -pp -threads 12 

done