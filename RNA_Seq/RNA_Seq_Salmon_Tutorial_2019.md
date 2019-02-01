---
title: "Differential Expression analysis using RNA-Seq data with DESeq2 (salmon)"
author: "Ian Dworkin"
date: "January 25th 2019"
output:
  html_document:
    keep_md: yes
    toc: yes
---


Modified from my NGS2016 Tutorial on Differential Expression analysis using RNA-Seq data and DESeq2


## Background
In this tutorial I will provide a basic overview of differential expression analysis for transcriptional profiling using RNA-Seq data. We will be using the [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) library in R. This approach utilizes a variant on the assumption of a negative binomially set of counts. This approach assumes that all you have going in are counts, that have not been normalized either for library size (or number of mapped reads), not for transcript length.

## Data provenance

The RNA-seq data we are using was generated from *Drosophila melanogaster*, where the developing wing and genital tissues (imaginal discs) were dissected out of larvae. We grew these flies at multiple temperatures (17C and 24C), where they tend to be bigger at lower temperatures. We also had "fed" and "starved" treatments (with starvation during development generally making organisms smaller). Importantly the wing tends to grow isometrically and is plastic with respect to nutrition and temperature. The genital discs much less so.

For each treatment combination, 3 independent biological samples (each sample consisting of ~30 imaginal discs) were produced. However several sequencing libraries failed, so the design is no longer balanced. 

So file names like

`samw_17C_gen_fed_R1_TAGCTT_L002`

Means that the genotype was Samarkand *white* (all of these samples are genetically identical), reared at 17 degrees C with high food, this would be replicate 1. L002 means this sample was run on lane 2 of the flow cell.

`samw_24C_wings_fed_R2_AGTCAA_L004` means the animals were reared at 24C, these were wings, also fed from lane 4. The 6 letter sequence is the barcode used for the sample during multiplexed sequencing.

In total, 20 samples were sequenced (100bp paired end using Illumina Tru-Seq chemistry)

## Checking the reads

All of the data for the class will be in.

`/2/scratch/Bio722_2019`

The data we are working with today is in

`/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset`

So navigate your way to there. take a look at the files using `ls`

Go to the `rawReads` directory

Question 1: How many files are in this directory? 
Question 2: What kind of files are these? 
Question 3: How many lines are there per distinct read in this file type?
Question 4: How many reads are there in each file (just do this for a few files, not all of them)
Question 5: If I had not unzipped each .fastq file, how would you have had to change your answer to question 4?


```bash
grep "^@HWI" samw_wings_starved_R3_GCCAAT_L004_R2_001.fastq | wc -l

#or

grep -c "^@HWI" samw_wings_starved_R3_GCCAAT_L004_R2_001.fastq 
```


## Quality control of raw reads

As always we start with doing some initial QC of the raw reads themselves. For this we will use `FastQC`. There are many other QC tools, but this is easy and fast. On info most bioinformatic software has been installed in `/usr/local`

In this case if you were to run fastqc to generate a report, you could do something like


```bash
fastqc YourFile.fastq -o YourOutputDirect/
```

Where you can specify which .fastq (or .fastq.gz) files and where you want to output them (`-o`). Note, please do all of your work in your own scratch folder! Not the one for the class. The fastqc has many options (as do most programs you are using) including specifying additional adapters and how to incorporate multiple .fastq files.

Question 6: Please run `fastqc` on a single .fastq file, but place the output in your own directory 

Question 7: While I don't want you to run it, how would run fastqc on each .fastq file sequentially?


```bash
fastqc *.fastq -o YourOutputDirect/
```

If you do this you may wish to specify the `-t n` for multi threads, with n being the number of processors.

We will take a look at one of the reports `fastqc` generated in class. 


## Trimming reads

Generally we would like to remove adapter sequences, and generally we want to consider removing low quality parts of reads as they are more likely to have incorrect sequences. There are many tools for this (scythe, cutadapt, BBMap), but we will use Trimmomatic, which is a fairly standard tool. As time is a bit short, we will not run this on full sized files. However here is an example script for one pair of files (a paired end run). Please do not run this!


```bash
dir=/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/rawReads
name=samw_wings_fed_R1_TGACCA_L002
dir_out=/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/trimmedReads

java -jar /usr/local/trimmomatic/trimmomatic-0.36.jar PE -threads 16 \
  ${dir}/${name}_R1_001.fastq ${dir}/${name}_R2_001.fastq \
  ${dir_out}/${name}_1_pe.fastq ${dir_out}/${name}_1_se.fastq \
  ${dir_out}/${name}_2_pe.fastq ${dir_out}/${name}_2_se.fastq \
  ILLUMINACLIP:/2/scratch/Bio722_2019/ID/BBMap_adapters.fa:2:30:10 LEADING:3 TRAILING:3 MAXINFO:20:0.2 MINLEN:36
```

So what have I done here? First I generated a few shell variables for input and output directories. I have also generated a shell variable for the name of the file pairs.

I am then running trimmomatic (which requires java) with a few options which I will explain a bit in class, but the manual is very useful so check that out as well.

We can run this one on the subsetted files. Again please set output directories into your own home directory or your own scratch directory.



```bash
dir=/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/rawReads/subset
name=subset_samw_17C_gen_fed_R1_TAGCTT_L002
dir_out=/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/trimmedReads

java -jar /usr/local/trimmomatic/trimmomatic-0.36.jar PE -threads 12 \
  ${dir}/${name}_R1_001.fastq ${dir}/${name}_R2_001.fastq \
  ${dir_out}/${name}_1_pe.fastq ${dir_out}/${name}_1_se.fastq \
  ${dir_out}/${name}_2_pe.fastq ${dir_out}/${name}_2_se.fastq \
  ILLUMINACLIP:/2/scratch/Bio722_2019/ID/BBMap_adapters.fa:2:30:10 LEADING:3 TRAILING:3 MAXINFO:20:0.2 MINLEN:36
```

The two most important things I want to go over today are the choice of adapter file and the parameters for MAXINFO. We will discuss these in class. 

Obviously we don't want to run the script each time so instead, I recommend making a for loop and running them all (consider submitting it as a shell script). It could look something like (PLEASE DON'T RUN THIS DURING THE TUTORIAL)


```bash
raw_dir=/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/rawReads
trim_dir=/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/trimmedReads


files=(${raw_dir}/*_R1_001.fastq)

for file in ${files[@]}
do
  name=${file}
  base=`basename ${name} _R1_001.fastq`
  java -jar /usr/local/trimmomatic/trimmomatic-0.36.jar PE -threads 16    ${raw_dir}/${base}_R1_001.fastq ${raw_dir}/${base}_R2_001.fastq  ${trim_dir}/${base}_R1_PE.fastq ${trim_dir}/${base}_R1_SE.fastq  ${trim_dir}/${base}_R2_PE.fastq ${trim_dir}/${base}_R2_SE.fastq   ILLUMINACLIP:/2/scratch/Bio722_2019/ID/BBMap_adapters.fa:2:30:10 LEADING:3 TRAILING:3 MAXINFO:20:0.2 MINLEN:36
 done
```

I will go through how this works, in particular `${files[@]}` in class.

All of the trimmed reads I generated are in 
`/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/trimmedReads`
Don't forget to perform fastqc (or whatever qc) on the trimmed files to confirm that you have removed adapteds as well as low quality sequence!

As mentioned in class, the trimmomatic adapter files are a bit minimal. I decided to use the adapter file from BBMap, which has many more adapters (including all of the index adapters in this data set). You should consider carefully which adapters you need to include. The fastqc report helps with this!


## How counts were generated
See Rob Patro's tutorial on using Salmon [here](http://angus.readthedocs.io/en/2016/rob_quant/tut.html). or [here](https://combine-lab.github.io/salmon/).


### The commands used for salmon for this data set
In case you want to try this yourself at a later date. DO NOT RUN this now.

First I downloaded the Drosophila transcriptome from flybase (in the drosophilaReference folder). 

salmon requires the generation of the index for the transcriptome (this only has to be done once per transcriptome). I used the commands

 **Don't re-run the index right now**

```bash
salmon index -t dmel-all-transcript-r6.25.fasta -i ./salmon_index/dmel-all-transcript_r6.25_index
```

With `-t` being the name of the transcriptome file being used, and `-i` specifying the location of the index being generated.

Once the index was generated I could then generate counts, using the trimmed paired end files. Here is an example of doing it for a set of paired end reads.


```bash
index_dir=/2/scratch/Bio722_2019/ID/drosophilaReference/salmon_index/dmel-all-transcript_r6.25_index

sample_dir=/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/trimmedReads

sample_name=samw_wings_starved_R3_GCCAAT_L004

out_dir=/2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/salmon_counts

salmon quant -i ${index_dir} -l A \
  -1 ${sample_dir}/${sample_name}_R1_PE.fastq \
  -2 ${sample_dir}/${sample_name}_R2_PE.fastq \
  -p 8 --validateMappings --rangeFactorizationBins 4 \
  --seqBias --gcBias \
  -o ${out_dir}/${sample_name}_quant
```

Work on writing your own for loop (with an array) to do the whole set. I recommend writing it as its own shell script.

I have a number of optional flags that I have set. These improve the mapping and quantification for both gcBias and other sequencing biases. These slow down the quantification a bit, but it is still typically less than 10 minutes for these samples (~30 million read pairs per sample)

Question 8: Trying running salmon on one sample (set of reads pairs) that differs from the ones in this example. As always, please output it in your own home directory.


### Getting the full set of counts we are going to use.
Counts from Salmon are found on info
`2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/salmon_counts`

I suggest using scp (with `-r`) to copy these to your local computer. Something like. There may be a better way, but I first generally copy what I am going to scp over to my regular folder

So on info (need to be logged into info11*, with * being the one you work on)


```bash
cd /2/scratch/Bio722_2019/ID/drosophilaDiscsGrowthSubset/

cp -r salmon_counts ~
```

Then on your laptop/local machine

```bash
cd YourWorkingDirHere
scp -r yourinfo@info.mcmaster.ca:~salmon_counts .
```

Which should copy the files to your local machine


## One other file we are going to need.

Salmon is mapping to transcripts, which is great, but we want to get gene level counts. As we will see the R library `tximport` has some utilities to help with this. However one thing we need is a file (normally named `txp_to_gene.tsv`) which is a tab delimited file with two columns (no header). The first column is going to be the name of the transcript, and the second column will be gene identifier. So the first few lines may look like this:



```bash
FBtr0005088	FBgn0260439
FBtr0006151	FBgn0000056
FBtr0070000	FBgn0031081
FBtr0070001	FBgn0052826
FBtr0070002	FBgn0031085
FBtr0070003	FBgn0062565
```

How do we generate this? In the R vignette for tximport `vignette('tximport')` it provides exaple code for this. You can potentially also get pretty close to a file like this for many model organisms. For *Drosophila* you can download a similar file from flybase that file link is here

ftp://ftp.flybase.org/releases/FB2018_06/precomputed_files/genes/fbgn_fbtr_fbpp_fb_2018_06.tsv.gz

I have already downloaded this onto info in `/2/scratch/Bio722_2019/ID/drosophilaReference
`. Go check it out.


```bash
head fbgn_fbtr_fbpp_fb_2018_06.tsv
tail fbgn_fbtr_fbpp_fb_2018_06.tsv
```

As you will see it has a few rows at the top and bottom of the file that you will want to remove, and we also want to remove the third column (FlyBase_FBpp are the protein identifiers). How would you do these to generate the file in the right format? You can use sed (or awk), but an even easier approach is


```bash
tail -n +7 fbgn_fbtr_fbpp_fb_2018_06.tsv | head -n -3 | cut -f 1,2 > temp.tsv
paste temp.tsv temp.tsv | cut -f 2,3 > txp_to_gene.tsv
rm temp.tsv

wc -l txp_to_gene.tsv
```

With the `tail -n +7` I am streaming all but the first 7 lines of the file which I pipe to `head -n -3` which streams all but the final 3 lines, but which is piped to `cut -f 1,2` which extracts just the first two columns. I created a temp file for us to take a look at to confirm it is what we expect. The only problem with this file is that the columns are in the reverse order.

So I use `paste` to make duplicate copy, and can extract the 2nd and 3rd column with cut

Or just use awk, which is a nice powerful way to do text processing.


```bash
awk '/^FBgn/ {print $2,"\t" $1}' fbgn_fbtr_fbpp_fb_2018_06.tsv > txp_to_gene_alt.tsv

wc -l txp_to_gene_alt.tsv
head txp_to_gene_alt.tsv
tail txp_to_gene_alt.tsv
```

For either of these files we may want to make sure there are not duplicate rows (i.e. cases where gene and transcript are the same but multiple proteins could be produced resulting in duplicate rows for what we care about)


```bash
sort -k2,2 -k1,1 txp_to_gene.tsv | uniq -d | less
```

Now we will move onto `R`
