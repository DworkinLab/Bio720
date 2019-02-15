# STAR for loop

#First set up some shell variables that we will want for directories, and files that do not change.

genome=/2/scratch/Bio722_2019/ID/drosophilaReference/star_index
trim_dir=/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/trimmedReads
mapped_dir=/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/star_mapped

#echo ${trim_dir}
#ls  ${trim_dir}

#echo ${genome}
#ls ${genome}

#echo ${mapped_dir}
#ls ${mapped_dir}



files=(${trim_dir}/*_R1_PE.fastq)
#echo ${files[@]}
#echo ${files[0]}

# Check that everything is working.
for file in ${files[@]}
do
    name=${file}
      echo ${name}
      base=`basename ${name} _R1_PE.fastq`
      echo ${base}
      echo ${base}_R1_PE.fastq
      echo ${base}_R2_PE.fastq
      echo ${name}
      echo ${trim_dir}/${base}_R1_PE.fastq
      echo ${trim_dir}/${base}_R2_PE.fastq
    /usr/local/STAR_2.6.1/bin/Linux_x86_64/STAR --version
done

# Note that (I assume for speed) you can keep the reference genome loaded
# into memory with STAR. So we want to do this before we run the loop,
# and then turn it off after we finish mapping everything.

#  load genome for shared use.
/usr/local/STAR_2.6.1/bin/Linux_x86_64/STAR --genomeLoad LoadAndExit \
  --genomeDir ${genome}


for file in ${files[@]}
do
  name=${file}
  base=`basename ${name} _R1_PE.fastq`
/usr/local/STAR_2.6.1/bin/Linux_x86_64/STAR --runThreadN 16 \
    --quantMode TranscriptomeSAM GeneCounts \
    --genomeDir ${genome} \
    --readFilesIn ${trim_dir}/${base}_R1_PE.fastq \
                  ${trim_dir}/${base}_R2_PE.fastq \
    --outFileNamePrefix ${mapped_dir}/${base} \
    --outSAMtype BAM SortedByCoordinate
done


# remove loaded genome

/usr/local/STAR_2.6.1/bin/Linux_x86_64/STAR --genomeLoad Remove \
  --genomeDir ${genome}



## Troubleshooting. It seems to not be finding the files, but I can not see what it is missing. Trying it with first pair of reads.
name=${files[0]}
echo ${name}
base=`basename ${name} _R1_PE.fastq`
echo ${base}
echo ${trim_dir}/${base}_R1_PE.fastq


## Could not diagnose the issue, let me try hard coding all

/usr/local/STAR_2.6.1/bin/Linux_x86_64/STAR --runThreadN 16 \
    --quantMode TranscriptomeSAM GeneCounts \
    --genomeDir /2/scratch/Bio722_2019/ID/drosophilaReference/star_index \
    --readFilesIn /2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/trimmedReads/samw_wings_starved_R1_CAGATC_L003_R1_PE.fastq \
    /2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/trimmedReads/samw_wings_starved_R1_CAGATC_L003_R2_PE.fastq \
    --outFileNamePrefix /2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/star_mapped/trial_samw_wings_starved_R1_CAGATC_L003_R2_PE.fastq \
    --outSAMtype BAM SortedByCoordinate

# just change genome dir - worked
# add trim_dir - worked
# add mapped_dir - worked
genome=/2/scratch/Bio722_2019/ID/drosophilaReference/star_index
trim_dir=/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/trimmedReads
raw_dir=/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/rawReads
mapped_dir=/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/star_mapped


/usr/local/STAR_2.6.1/bin/Linux_x86_64/STAR --runThreadN 16 \
    --quantMode TranscriptomeSAM GeneCounts \
    --genomeDir ${genome} \
    --readFilesIn ${trim_dir}/samw_wings_starved_R1_CAGATC_L003_R1_PE.fastq \
                  ${trim_dir}/samw_wings_starved_R1_CAGATC_L003_R2_PE.fastq \
    --outFileNamePrefix ${mapped_dir}/trial_samw_wings_starved_R1_CAGATC_L003_QC \
    --outSAMtype BAM SortedByCoordinate

# also for QC use the raw .fastq files, not the trimmed files. Like this.
# It only wants one primary read.
/usr/local/STAR_2.6.1/bin/Linux_x86_64/STAR --runThreadN 16 \
    --quantMode TranscriptomeSAM GeneCounts \
    --genomeDir ${genome} \
    --readFilesIn ${raw_dir}/samw_wings_starved_R3_GCCAAT_L004_R1_001.fastq \
                  ${raw_dir}/samw_wings_starved_R3_GCCAAT_L004_R2_001.fastq \
    --outFileNamePrefix ${mapped_dir}/trial_samw_wings_starved_R1_CAGATC_L003_QC \
    --outSAMtype BAM SortedByCoordinate
