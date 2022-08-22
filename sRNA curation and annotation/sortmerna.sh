#!/bin/bash 

for i in `ls -1 *_trim.fq.gz | sed 's/_trim.fq.gz//'` 

do

echo "Uncompressing FASTQ data of ${i}"
gunzip ${i}_trim.fq.gz

echo "Running SortMeRNA for ${i}"
/path/to/sortmerna --ref all_rRNA_databases --reads "${i}_trim.fq" --aligned "${i}_rRNA" --other "${i}_no_rRNA" --log -a 32 -v --fastx

done

##--- END HERE --- 
