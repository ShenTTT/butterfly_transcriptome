#!/bin/bash 

for i in `ls -1 *_1P.fq.gz | sed 's/_1P.fq.gz//'` 

do 

echo "Perform hisat2 for ${i}"

/path/to/hisat2 -p 32 --summary-file ${i}.txt -x genome -1 ${i}_1P.fq.gz -2 ${i}_2P.fq.gz -S ${i}_map.sam

echo "convert sam to bam for ${i}"
/path/to/samtools view -@ 24 -bS ${i}_map.sam > ${i}_map.bam

echo "sort bam file of ${i}"
/path/to/samtools sort ${i}_map.bam -o ${i}_final.bam -@ 24

echo "index ${i}"
/path/to/samtools index ${i}_final.bam

done
