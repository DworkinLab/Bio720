#!/bin/bash -login

#PBS -l walltime=02:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=4gb
#PBS -M MyEmail@server.edu
#PBS -N sd_BSA_2014_trimmomatic
#PBS -j oe
#PBS -m abe
#PBS -t 0-7

module load java
module load Trimmomatic
cd $PBS_O_WORKDIR
dir=/mnt/scratch/idworkin/sd_BSA/

#Calculate File names

#make an array for each file in the directory "dir" that ends in _R1.fastq
files=(${dir}*_R1.fastq.gz)

# i.e. use the echo command to see how each element of the array
#echo ${files[0]}
#echo ${files[1]}

# pulls out one element at a time of the array based on PBS_ARRAYID (where PBS_ARRAYID is a special variable for the integers in PBS -t above.
file=${files[${PBS_ARRAYID}]}

# creates variable "name" for the current file based on "file"
# basename is a program that strips folder and suffix information from filenames.
name=`basename $file _R1.fastq.gz`

# regex to get first part of file name (excludes the extension)
bname=${name%.*}

#check to make sure that it produces the correct command
#echo "wc -l ${dir}${name} > $bname.out"

# run it!
#wc -l ${dir}${name} > $bname.out
java -jar $TRIM/trimmomatic PE ${dir}${name}_R1.fastq.gz ${dir}${name}_R2.fastq.gz \
 ${dir}trimmed/${name}_R1_pe ${dir}trimmed/${name}_R1_se ${dir}trimmed/${name}_R2_pe \
 ${dir}trimmed/${name}_R2_se \
 ILLUMINACLIP:/mnt/research/common-data/Bio/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
 LEADING:3 TRAILING:3 MAXINFO:40:0.3 MINLEN:40

qstat -f $PBS_JOBID
