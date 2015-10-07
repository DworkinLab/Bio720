#!/bin/bash

files=”Heart
Liver
Lung
Muscle
Kidney
Ovary”

for file in $files
do
	java -jar /work/ben/Trimmomatic-0.32/trimmomatic-0.32.jar PE -phred33 -trimlog ${file}_log.txt AO248_${file}_R1.fastq.gz AO248_${file}_R2.fastq.gz AO248_${file}_R1_trim_paired.fastq.gz AO248_${file}_R1_trim_single.fastq.gz AO248_${file}_R2_trim_paired.fastq.gz AO248_${file}_R2_trim_single.fastq.gz ILLUMINACLIP:/work/ben/Trimmomatic-0.32/adapters/trimmomatic_adapters.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36
done
