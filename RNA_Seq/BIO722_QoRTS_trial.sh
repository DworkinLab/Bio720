
# variables
mapped_dir=/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/star_mapped
qc_dir=/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/QC_reports/QC_STAR

bam_file=trial_samw_wings_starved_R1_CAGATC_L003_QCAligned.sortedByCoord.out.bam

base=`basename ${bam_file} .sortedByCoord.out.bam`

gtf_file=/2/scratch/Bio722_2019/ID/drosophilaReference/dmel-all-r6.25.gtf

genome_fasta=/2/scratch/Bio722_2019/ID/drosophilaReference/dmel-all-chromosome-r6.25.fasta

raw_dir=/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/rawReads
raw1=samw_wings_starved_R3_GCCAAT_L004_R1_001.fastq
raw2=samw_wings_starved_R3_GCCAAT_L004_R2_001.fastq

# command
mkdir ${qc_dir}/${base}
java -jar -Xmx8G /usr/local/qorts/QoRTs.jar QC \
                   --generatePlots \
                   --genomeFA ${genome_fasta} \
                   --rawfastq ${raw_dir}/${raw1},${raw_dir}/${raw1} \
                   ${mapped_dir}/${bam_file} \
                   ${gtf_file} \
                   ${qc_dir}/${base}/
