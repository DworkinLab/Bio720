# Mapping and counting transcripts using RNASeq, a basic introduction

This tutorial is based upon a tutorial I ran with Dr. Chris Chandler at the NGS2015 workshop. The original tutorial can be found [here](http://angus.readthedocs.org/en/2015/drosophila_rnaseq_counting_htseq.html).

## Background
In this tutorial, we’ll use some sample data from a project we did on flies (*Drosophila melanogaster*) to illustrate how you can use RNA-Seq data to look for differentially expressed genes. Here’s some brief background on the project: we’re trying to understand how different wild-type genetic backgrounds can influence the phenotypic effects of mutations, using the developing fly wing as our model system. We have several mutations that disrupt wing development, and we’ve backcrossed them into the genetic backgrounds of two different wild-type fly strains (Samarkand/SAM and Oregon-R/ORE). For this tutorial, we’ve taken a subset of the data–we’ll look for expression differences in developing wing tissues between non-mutant and mutant (sd[E3]) flies in each of the two genetic backgrounds, and in flies with a “hybrid” genetic background (i.e., crosses between SAM/ORE flies, again both with and without the mutation). To make things run a little bit faster, we’ve included only 250K reads for each sample.


## Connecting to a remote machine using `ssh`

At this point you should be able to do this in your sleep...

## The files we are working with.

The files we want should be in
```
/net/infofile2/2/scratch/Bio720_ID/drosophilaReads
 ```


## Making a useful directory structure

For today, we can actually house all of these smaller files in your home directory. However, normally it would be best to keep them in scratch for full sequencing runs.

Let's make some directories starting in your home directory

```
cd

mkdir drosophila_rnaseq
cd drosophila_rnaseq
mkdir {reads,trimmed,mapped,counts,references}
```

Now copy the reads over to the directory `trimmed`

```
cp /net/infofile2/2/scratch/Bio720_ID/drosophilaReads/*.fastq  ~/drosophila_rnaseq/trimmed
```

## Tools we are using
 For today we will take a basic approach and will map reads using one of the well known splice aware aligners [tophat](https://ccb.jhu.edu/software/tophat/index.shtml), and use [HTseq](http://www-huber.embl.de/HTSeq/doc/overview.html) to get counts. As we discussed in class there are many different approaches to mapping reads to genomes/transcriptomes and to generating counts. Thankfully Brian has already installed the tools we are using today on his cluster. I do want to say that I chose these for ease of use (since you are just learning), and other tools may be much better suited to your needs.

## Preparing and indexing the reference.
 Next, we need to get our reference genome. This is another area where you want to be careful and pay attention to version numbers–the public data from genome projects are often updated, and gene ID’s, coordinates, etc., can sometimes change. At the very least, you need to pay attention to exactly which version you’re working with so you can be consistent throughout all your analyses.

In this case, we’ll download the reference Drosophila genome (`.fa`) and the `.gtf` annotation file (which has the ID’s and coordinates of known transcripts, etc.) from ensembl. We’ll put it in its own directory so we keep our files organized. Then we’ll prepare the genomes for use with our software tools by indexing them:

NOTE.. THese steps take a few minutes, so you may wish to skip them and use the pre-computed ones for now.

```
cd ~/drosophila_rnaseq/references

curl -O ftp://ftp.ensembl.org/pub/release-75/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP5.75.dna.toplevel.fa.gz
gunzip Drosophila_melanogaster.BDGP5.75.dna.toplevel.fa.gz
curl -O ftp://ftp.ensembl.org/pub/release-75/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP5.75.gtf.gz
gunzip Drosophila_melanogaster.BDGP5.75.gtf.gz

bowtie-build Drosophila_melanogaster.BDGP5.75.dna.toplevel.fa Drosophila_melanogaster5_75
samtools faidx Drosophila_melanogaster.BDGP5.75.dna.toplevel.fa
```

Running this block of code may take a bit of time (especially the bowtie-build command). If you don’t want to wait for it, I have included some pre-computed reference files for you.

## Mapping reads
Now we’re ready for the first step: mapping our RNA-seq reads to the genome. We will use tophat+bowtie1, which together are a splicing-aware read aligner.

###QC and trimmming
Don’t forget that with your reads, you’ll probably want to take care of the usual QC steps before you actually begin your mapping. The cleaned_reads directory contains reads that have already been filtered (adapter removal, etc.) and lightly trimmed. Just as a reminder, for a pair of samples you could run `trimmomatic` as follows: (DON'T RUN)
```
#java -jar /usr/local/trimmomatic/trimmomatic-0.33.jar PE -threads 8 \
#SAM_w_r1_ACTTGA_L003_R1_001.fastq \
#SAM_w_r1_ACTTGA_L003_R2_001.fastq \
#../trimmed/SAM_wt_rep1_PE_1.fastq ../trimmed/SAM_wt_rep1_SE_1.fastq \
#../trimmed/SAM_wt_rep1_PE_2.fastq ../trimmed/SAM_wt_rep1_SE_2.fastq \
#ILLUMINACLIP:/usr/local/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10
```

As it turns out, for many of the read pairs (since we are only looking at 250K reads for each sample), there are very few cases where we have singletons left after trimming (i.e. only one of the two reads from a cluster survived trimming). Thus for the moment we are ignoring the residual singletons. It is up to you for own work, whether it is worth utilizing the reads from the singletons. It can a little bit of pain to the scripts, but most tools accept them and can map them along with the surviving paired end reads.

###Back to mapping
Since we have a lot of files to map, it would take a long time to re-write the mapping commands for each one. And with so many parameters, we might make a mistake or typo. It’s usually safer to use a simple shell script with shell variables to be sure that we do the exact same thing to each file. Using well-named shell variables also makes our code a little bit more readable:


We will create shell variables to store the location of our reference genome and annotation file. Note that we are leaving off the .fa from the reference genome file name, because some of the later commands will require just the base of the file name.

```
reference=~/drosophila_rnaseq/references/Drosophila_melanogaster.BDGP5.75.dna.toplevel
annotation=~/drosophila_rnaseq/references/Drosophila_melanogaster.BDGP5.75.gtf
```

We can check quickly with an echo command
```
echo $reference
echo $annotation
```

Now Make sure we are in the right directory. We will store all of our mapped reads in the `mapped` sub-folder in our project folder

```
cd ~/drosophila_rnaseq
```

Now we Create an array to hold the names of all our samples. Later, we can then cycle through each sample using a simple for loop. We are actually only doing a subset of all of the files. However, I will provide all of the counts for the full data set at the end, so we can play with it.

```
samples[1]=ORE_wt_rep1_PE
samples[2]=ORE_wt_rep2_PE
samples[3]=ORE_sdE3_rep1_PE
samples[4]=ORE_sdE3_rep2_PE
samples[5]=SAM_wt_rep1_PE
samples[6]=SAM_sdE3_rep1_PE
samples[7]=SAM_sdE3_rep2_PE
samples[8]=HYB_sdE3_rep1_PE
samples[9]=HYB_sdE3_rep2_PE
```

We can check the array variable we have just created

## Now we can actually do the mapping

```
for i in 1 2 3 4 5 6 7 8 9
do
    sample=${samples[${i}]}
    echo ${sample}
    #Map the reads
    /usr/local/tophat-2.0.8/tophat -p 4 -G ${annotation} -o ~/drosophila_rnaseq/mapped/${sample} ${reference} \
     ~/drosophila_rnaseq/trimmed/${sample}_1.fastq \
     ~/drosophila_rnaseq/trimmed/${sample}_2.fastq \

    #Count the number of reads mapping to each feature using HTSeq
    htseq-count --format=bam --stranded=no --order=pos ~/drosophila_rnaseq/mapped/${sample}/accepted_hits.bam ${annotation} > ~/drosophila_rnaseq/counts/${sample}_htseq_counts.txt
done
```

## The counts
If you navigate to the `./counts` sub-directory (inside the project folder), you will see we have now generated one count file for each of the pairs of samples. Let's take a look at these
