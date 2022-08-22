#!/bin/bash 

for i in `ls -1 *_trim.fq.gz | sed 's/_trim.fq.gz//'` 

do

echo "Uncompressing FASTQ data of ${i}"
gunzip ${i}_trim.fq.gz

echo "Running SortMeRNA for ${i}"
/home/nus/Desktop/NGS/sortmerna-2.1b/sortmerna --ref /home/nus/Desktop/NGS/sortmerna-2.1b/rRNA_databases/silva-bac-16s-id90.fasta,/home/nus/Desktop/NGS/sortmerna-2.1b/index/silva-bac-16s-db:/home/nus/Desktop/NGS/sortmerna-2.1b/rRNA_databases/silva-bac-23s-id98.fasta,/home/nus/Desktop/NGS/sortmerna-2.1b/index/silva-bac-23s-db:/home/nus/Desktop/NGS/sortmerna-2.1b/rRNA_databases/silva-arc-16s-id95.fasta,/home/nus/Desktop/NGS/sortmerna-2.1b/index/silva-arc-16s-db:/home/nus/Desktop/NGS/sortmerna-2.1b/rRNA_databases/silva-arc-23s-id98.fasta,/home/nus/Desktop/NGS/sortmerna-2.1b/index/silva-arc-23s-db:/home/nus/Desktop/NGS/sortmerna-2.1b/rRNA_databases/silva-euk-18s-id95.fasta,/home/nus/Desktop/NGS/sortmerna-2.1b/index/silva-euk-18s-db:/home/nus/Desktop/NGS/sortmerna-2.1b/rRNA_databases/silva-euk-28s-id98.fasta,/home/nus/Desktop/NGS/sortmerna-2.1b/index/silva-euk-28s:/home/nus/Desktop/NGS/sortmerna-2.1b/rRNA_databases/rfam-5s-database-id98.fasta,/home/nus/Desktop/NGS/sortmerna-2.1b/index/rfam-5s-db:/home/nus/Desktop/NGS/sortmerna-2.1b/rRNA_databases/rfam-5.8s-database-id98.fasta,/home/nus/Desktop/NGS/sortmerna-2.1b/index/rfam-5.8s-db --reads "${i}_trim.fq" --aligned "${i}_rRNA" --other "${i}_clean" --log -a 32 -v --fastx

echo "Doing cleanup for ${i}"
gzip ${i}_trim.fq ${i}_clean.fq ${i}_rRNA.fq


done

##--- END HERE --- 
