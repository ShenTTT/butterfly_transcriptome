#!/bin/bash 

for i in `ls -1 *_1P.fq.gz | sed 's/_1P.fq.gz//'` 

do 

echo "Run salmon for ${i}"

salmon quant -i /path/to/gentrome_index -p 36 --libType A -g /path/to/genome.gtf -1 ${i}_1P.fq.gz -2 ${i}_2P.fq.gz --validateMappings --seqBias --gcBias -o /output/to/${i}_quant

done
