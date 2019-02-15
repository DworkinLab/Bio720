# STAR

First generate the genome index

/usr/local/STAR_2.6.1/bin/Linux_x86_64/STAR --runThreadN 12 --runMode genomeGenerate \
  --genomeDir /2/scratch/Bio722_2019/ID/drosophilaReference/star_index \
  --genomeFastaFiles /2/scratch/Bio722_2019/ID/drosophilaReference/dmel-all-chromosome-r6.25.fasta \
  --sjdbGTFfile /2/scratch/Bio722_2019/ID/drosophilaReference/dmel-all-r6.25.gtf \
  --sjdbOverhang 100 --genomeSAindexNbases 12


/usr/local/STAR_2.6.1/bin/Linux_x86_64/STAR --runThreadN 12 \
    --genomeDir /2/scratch/Bio722_2019/ID/drosophilaReference/star_index \
    --readFilesIn /2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/trimmedReads/subset_samw_17C_gen_fed_R1_TAGCTT_L002_1_pe.fastq \
    /2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/trimmedReads/subset_samw_17C_gen_fed_R1_TAGCTT_L002_2_pe.fastq \
    --outFileNamePrefix /2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/star_mapped/subset_samw_17C_gen_fed_R1 \
    --outSAMtype BAM SortedByCoordinate
