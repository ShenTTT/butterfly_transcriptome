#!/bin/bash 

for i in `ls -1 *.fq.gz | sed 's/.fq.gz//'` 

do 

echo "Run trimmomatic for ${i}"

java -jar /path/to/trimmomatic.jar SE -threads 24 $i\.fq.gz $i\_trim.fq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10

done

##--- END HERE --- 