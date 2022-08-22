#!/bin/bash 

for i in `ls -1 *_trim.fq | sed 's/_trim.fq//'` 
do
echo "mapping ${i}"
/home/nus/Desktop/NGS/STAR-2.7.8a/bin/Linux_x86_64/STAR --runThreadN 32 --alignIntronMax 1 --outFilterMultimapNmax 100000 --outFilterMismatchNmax 3 --readFilesIn ${i}_trim.fq --genomeDir /media/nus/LaCie/BANY_v1.2/index/STAR_index --outFileNamePrefix ./${i}_

done
