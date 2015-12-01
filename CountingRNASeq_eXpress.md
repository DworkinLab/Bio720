# Mapping and counting transcripts using RNASeq data, with eXpress

In the first tutorial we did counted mapped reads to the genome (using a splice aware mapping tool, tophat) using HTSeq to perform the counting. One alternative approach, especially for non-model systems is to use contigs from a transcriptome to map the reads onto. In particular if you have RNA-Seq data from a non-model organism too divergent from (see papers...) a good reference genome (with good gene models) to map to may not be avilable. In addition some tools only accept mappings from a transcriptome.

Today instead of running this live, for the sake of time, we will just go over the tutorial and you can run this on your own outside of class.


## Background
We will be using the same data set from *Drosophila melanogaster* wing imaginal discs as we did in the previous tutorial. However this time we will be using a different piece of software for counting, namely [eXpress](http://bio.math.berkeley.edu/eXpress/). One great thing about `eXpress` is that it is, as with many great pieces of UNIX software, a **streaming tool** (see page 41-42 in your BDS text for an explanation of streaming). THis means it is fairly light on memory (unlike many bioinformatics tools), and in principle could be run on your UNIX laptops (linux, mac OS X, etc). However, it will be much faster to do it on a cluster where we have dedicated resources.In addition, if disk space is an image, `eXpress` accepts input on a stream as well. This means you can stream directly from the alignment tool (they suggest `bowtie` or `bowtie2`) into `eXpress` and only keep the output file.

If you plan on using `eXpress` for your own work, I highly recommend working through their [tutorial](http://bio.math.berkeley.edu/eXpress/tutorial.html) as well. In addition there is a useful tutorial on using this as part of a larger pipeline from Harold Pimentel available [here](http://lmcb.wikispaces.com/eXpress+Walkthrough).


## The files we are working with.

The files we want should be in.
```bash
/net/infofile2/2/scratch/Bio720_ID/drosophilaReads
 ```


## Make a new directory for the mapped reads and counts.

Keep in mind that when we did our previous mappings (using tophat) to the genome. So we want a fresh directory to work with. However, since these are both part of the same project, it is worth keeping them organized together.

Navigate back to the `drosophila_rnaseq` folder you made for the previous tutorial. Now we are going to make several new sub-directories, one for the new mapped reads, and one for the counts.


```bash
cd
cd drosophila_rnaseq
mkdir {transcriptome_mapped,counts_eXpress}
```

Reads that we will use should already be in the previously created `trimmed` folder


## Tools we are using

As discussed above we will be using `eXpress` for counting mapped reads. In addition to my comments above (and as always I recommend going through the tutorial and manual of those who wrote such tools) there are a few important things to keep in mind.
1. `eXpress` needs to keep track of multi-mapped reads. As such the authors have some recommended settings for the alignment process.  If you are using bowtie2 they suggest the following parameters during the call to the software.

```bash
-a -X 600 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 --no-discordant --no-mixed
```

2. As we discussed when we did the tutorial on using DESeq2, the output from `eXpress` may not always be in the same order. So make sure in your R script (or shell script) to sort by your gene names prior to any further analyses!

## Preparing and indexing the reference.

Again, as a reminder, we are not using the whole genome (and gene annotations) as a reference. Instead we are using the *Drosophila melanogaster* transcriptome. So we need to download a reference transcriptome (or make a *de novo* transcriptome assembly), and then prepare and index the reference.

NOTE.. THese steps takes quite a while, so you may wish to skip them and use the pre-computed ones for now.

```bash
cd ~/drosophila_rnaseq/references

curl -O ftp://ftp.flybase.net/releases/FB2015_05/dmel_r6.08/fasta/dmel-all-transcript-r6.08.fasta.gz
gunzip dmel-all-transcript-r6.08.fasta.gz
```

Now we can index this for the aligner (this is the step that can take 30+ minutes)

```bash
bowtie2-build --offrate 1 dmel-all-transcript-r6.08.fasta dmel-all-transcript-r6.08_bowtie2.index
```

Running this block of code may take a bit of time (especially the bowtie2-build command). If you don’t want to wait for it, I have included some pre-computed reference files for you. These can be found in...

## Mapping reads
Now we’re ready for mapping the reads to the transcriptome.

###QC and trimmming
Don’t forget that with your reads, you’ll probably want to take care of the usual QC steps before you actually begin your mapping. The cleaned_reads directory contains reads that have already been filtered (adapter removal, etc.) and lightly trimmed. Just as a reminder, for a pair of samples you could run `trimmomatic` as follows: (DON'T RUN). In this case I did this for the one sample (paired) that we are interested in.
```bash
java -jar /usr/local/trimmomatic/trimmomatic-0.33.jar PE -threads 8 \
ORE_w_r1_ATCACG_L001_R1_001.fastq \
ORE_w_r1_ATCACG_L001_R2_001.fastq \
./trimmed/ORE_wt_rep1_PE_1.fastq ./trimmed/ORE_wt_rep1_SE_1.fastq \
./trimmed/ORE_wt_rep1_PE_2.fastq ./trimmed/ORE_wt_rep1_SE_2.fastq \
ILLUMINACLIP:/usr/local/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10
```

After a few minutes (5ish with 8 threads going). you should see output like:

```bash
## Using PrefixPair: 'TACACTCTTTCCCTACACGACGCTCTTCCGATCT' and 'GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT'
## ILLUMINACLIP: Using 1 prefix pairs, 0 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
##Quality encoding detected as phred33
##Input Read Pairs: 31185254 Both Surviving: 30954796 (99.26%) Forward Only Surviving: 225488 (0.72%) Reverse Only Surviving: 0 (0.00%) Dropped: 4970 (0.02%)
##TrimmomaticPE: Completed successfully
```

For the moment we are ignoring the residual singletons (SE). IN this case this is ~225K singleton reads (how did I figure this out). In comparison for the remaining paired end reads there are 30954796 reads, so over 99% survived. This is in part because I mostly did adapter trimming (which is ok for this, why?). However, `eXpress` can handle the combination of paired and single end reads. It is up to you for own work, whether it is worth utilizing the reads from the singletons. It can a little bit of pain to the scripts, but most tools accept them and can map them along with the surviving paired end reads.

###Back to mapping
For this tutorial, I am only going to use one set of reads (paired end) from one sample. Please check out the previous tutorial for an example of how you can use a `for` loop to run though all of the samples.

We will create shell variables to store the location of our reference genome and annotation file. Note that we are leaving off the .fa from the reference genome file name, because some of the later commands will require just the base of the file name. I am creating some simple shell variables just to make running the commands clear.

Also note that for a couple of the variables, I have only used partial names. This is on purpose, and you will see why below.

```bash
ref_transcriptome=~/drosophila_rnaseq/references/dmel-all-transcript-r6.08
sample_PE=ORE_wt_rep1_PE
dir_in=~/drosophila_rnaseq/trimmed_full
dir_out=~/drosophila_rnaseq/transcriptome_mapped
```

We can check quickly with an echo command
```bash
echo ${ref_transcriptome}
```

You should try the rest yourself to make sure the variables we have created are the file or directory names that match appropriate locations!

We will store all of our mapped reads in the `transcriptome_mapped` sub-folder in our project folder

```
cd ~/drosophila_rnaseq
```

Now we Create an array to hold the names of all our samples. Later, we can then cycle through each sample using a simple for loop. We are actually only doing a subset of all of the files. However, I will provide all of the counts for the full data set at the end, so we can play with it.


We can check the array variable we have just created

Now we can actually do the mapping. The `-p 16` flag is how many threads we want to use. If we are all doing it at the same time, I suggest using 4 instead.

```bash
bowtie2 -a -X 1200 -p 16 --rdg 6,5 --rfg 6,5 --score-min L,-.6,-.4 --no-discordant --no-mixed \
   -x ${ref_transcriptome}.index \
   -1 ${dir_in}/${sample_PE}_1.fastq -2 ${dir_in}/${sample_PE}_2.fastq | samtools view -Sb - >  ${dir_out}/${sample_PE}_hits.bam
```
This will likely take a while to run.

Some of the flags we have set:
`-a` as we want it to report all alignments
`-X` setting the maximum fragment length to 1200bp
`-x` the bowtie2 index file we created
`--rdg` read gap penalty
`--rfg` reference gap penalty
`--score-min` sets options for valid alignments
`--no-discordant` excludes pairs that are not concordant (are not both contained in the same contig)
`--no-mixed` only includes matches for full pairs

These are the recommend `bowtie2` parameters as suggested by the authors of `eXpress`. These include limiting the length of indels and avoiding unlikely splice variants. More info about the individual parameters can be found [here](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)


## counting with `eXpress`
Once we have the .bam files, eXpress is actually pretty quick.

keep in mind that express does not allow you to name each output file, just output directories. So we have to do some moving around after we finish generating the counts. This is actually the main reason I tend not to use `eXpress`.


```bash
file_bam=~/drosophila_rnaseq/transcriptome_mapped/ORE_wt_rep1_PE
ref=~/drosophila_rnaseq/references/dmel-all-transcript-r6.08.fasta
output_count=~/drosophila_rnaseq/counts_eXpress

express -o ${output_count}/${file_bam}/ ${ref} ${file_bam}_hits.bam
```

`eXpress` does not allow renaming of the two output files `params.xprs` & `results.xprs`. So we rename them ourselves.

```
mv ${output_count}/${file_bam}/params.xprs ${output_count}/${file_bam}_params.xprs
mv ${output_count}/${file_bam}/results.xprs ${output_count}/${file_bam}_results.xprs
```
## The counts
If you navigate to the `~/drosophila_rnaseq/counts_eXpress` sub-directory (inside the project folder), you will see we have now generated one count file (results) for each of the pairs of samples. Let's take a look at these. All of these files should have the same number of rows, how can you check?  

 We can also ask how many "counts" we have. We should look at the output file first to figure out which column we want.
```bash
awk '{s+=$2} END{print s}' PickAFile
```
