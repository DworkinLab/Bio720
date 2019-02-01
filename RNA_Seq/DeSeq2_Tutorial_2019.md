---
title: "Differential Expression analysis using RNA-Seq data with DESeq2 (salmon)"
author: "Ian Dworkin"
date: "January 25th 2019"
output:
  html_document:
    keep_md: yes
    number_sections: yes
    toc: yes
---
 
 

 Modified from my NGS2016 Tutorial on Differential Expression analysis using RNA-Seq data and DESeq2


## Background
In this tutorial I will provide a basic overview of differential expression analysis for transcriptional profiling using RNA-Seq data. We will be using the [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) library in R. This approach utilizes a variant on the assumption of a negative binomially set of counts. This approach assumes that all you have going in are counts, that have not been normalized either for library size (or number of mapped reads), not for transcript length.

Instead of running these analyses on info, we'll run this locally on our own computers. Before you begin, you will need to download all of the count files we generated using Salmon. Rob Patro and the author of DESeq2 (Mike Love) have developed some nice import tools to get everything into `R` relatively efficiently. However, you will need to use a very recent version of `R` to use these functions.

## Install DeSeq2
Please install DeSeq2 in R or R studio on your laptop

## Data provenance

The RNA-seq data we are using was generated from Drosophila melanogaster, where the developing wing and genital tissues (imaginal discs) were dissected out of larvae. We grew these flies at multiple temperatures (17C and 24C), where they tend to be bigger at lower temperatures. We also had "fed" and "starved" treatments (with starvation during development generally making organisms smaller). Importantly the wing tends to grow isometrically and is plastic with respect to nutrition and temperature. The genital discs much less so.

For each treatment combination, 3 independent biological samples (each sample consisting of ~30 imaginal discs) were produced. However several sequencing libraries failed, so the design is no longer balanced. 

So file names like

`samw_17C_gen_fed_R1_TAGCTT_L002`

Means that the genotype was Samarkand *white* (all of these samples are genetically identical), reared at 17 degrees C with high food, this would be replicate 1. L002 means this sample was run on lane 2 of the flow cell.

`samw_24C_wings_fed_R2_AGTCAA_L004` means the animals were reared at 24C, these were wings, also fed from lane 4. The 6 letter sequence is the barcode used for the sample during multiplexed sequencing.

In total, 20 samples were sequenced (100bp paired end using Illumina Tru-Seq chemistry)


## How counts were generated
See Rob Patro's tutorial on using Salmon [here](http://angus.readthedocs.io/en/2016/rob_quant/tut.html). or [here](https://combine-lab.github.io/salmon/).


### The commands used for salmon for this data set
In case you want to try this yourself at a later date. DO NOT RUN This now.

First I downloaded the Drosophila transcriptome (in the drosophilaReference folder). 

salmon requires the generation of the index for the transcriptome (this only has to be done once per transcriptome). I used the commands

 **Don't re-run the index right now**

```bash
salmon index -t dmel-all-transcript-r6.25.fasta -i ./salmon_index/dmel-all-transcript_r6.25_index
```

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


## Get `R` loaded, and let's get started.

First we load in libraries in `R`. In addition to `DESeq2`, there are a few other libraries we will need.


```r
library(DESeq2)
library(tximport)
library(readr)
library("RColorBrewer")
library(gplots)
```

It is possible or even likely that you will get an error for some of these, as you have not yet installed the appropriate library. Some are from CRAN (where most R libraries are available), while others are part of bioconductor.

To install for base R (like `gplots`) you can use the:


```r
install.packages("gplots")
```

While for `tximport` you will need to do this via bioconductor:

Note only run this if tximport did not install!

```r
source("http://bioconductor.org/biocLite.R")
biocLite("tximport")
```

Depending on the implementation of R you are using (R on mac, R on windows, R studio), there may be some slight differences, so grab an instructor.

Set the working directory for the raw count data (you need to know where you put it). I will go over how I organize my projects to keep this simple.  Just like Unix, R has a current working directory. You can set the working directory using the `setwd()` function. Once you have unzipped the file Rob has provided, you need to navigate to that directory. You will want to be inside the quantification folder. For me this will look like this.

The folder that contains all of the sub-folders with the quantifications should be renamed "quants"


```r
#setwd("../data/salmon_counts")
# This will differ for you!!!
setwd("~/TeachingAndLectures/Bio722/Bio722_2019/salmon_counts/")
# Setting it up for the import (this is not the import itself)
quant_files <- file.path("quants", list.files("quants"), "quant.sf")

# Let's take a look at this variable we have created
quant_files
```

## Loading the count data into R

DESeq2 and other libraries often have helper functions for getting your count data in. In particular if you are using objects created from other tools that the same authors generated. However, if you are going to make your own pipeline, it is important to know how to write some simple R code to be able to get your data in, and to parse it so that it is in the format you need. I will (if we have time) go through a more typical example where there is no helper functions (so you write it yourself). However, we will use the ones available from `tximport`.


### Getting the meta-data set

Normally I extract these direct from the file names (which I could here given how I have the names written). However, more generally (but more laborious) you can do it directly in `R` like I have below, or you could make a second data frame containing the meta data.


```r
# Names of samples.
samples <- c("samw_17C_gen_fed_R1", 
             "samw_17C_gen_fed_R2", 
             "samw_17C_gen_fed_R3",
             "samw_17C_gen_starved_R1", 
             "samw_17C_gen_starved_R2", 
             "samw_17C_wings_fed_R1", 
             "samw_17C_wings_fed_R2",
             "samw_17C_wings_fed_R3", 
             "samw_17C_wings_starved_R1",
             "samw_17C_wings_starved_R2",
             "samw_17C_wings_starved_R3", 
             "samw_24C_gen_fed_R1",
             "samw_24C_gen_fed_R3",
             "samw_24C_gen_starved_R3",
             "samw_24C_wings_fed_R1",
             "samw_24C_wings_fed_R2",
             "samw_24C_wings_fed_R3",
             "samw_24C_wings_starved_R1",
             "samw_24C_wings_starved_R2",
             "samw_24C_wings_starved_R3")

names(quant_files) <- samples
```

We also need to load in the table of gene and transcript names (with correspondence).  We can first take a look at this file (`txp_to_gene.tsv`) in your terminal. How do you do this?

To import it into `R` we are going to use the most basic and general import tool in base `R`. `read.table()`. There are many options (and I often use read.csv). The data.table library has many useful approaches to this problem (FYI.)



```r
tx2gene <- read.table("txp_to_gene.tsv", col.names=c("TXNAME", "GENEID"))

head(tx2gene)
dim(tx2gene)
```
Drosophila genomics and datasets have their own unique gene and transcript identifiers (as well as more universal ones). In this case the column TXNAME are the identifiers for the transcripts (FBtr) and the gene identifier is in the GENEID column (FBgn). The names for those columns are from when we generated above for the counts using salmon. More information about Drosophila genes, transcripts, etc can be found at [flybase](http://flybase.org/). 

It is probably (before we start) worth asking about how many genes and how many transcripts we have in our list. How might you do this in `R`?


```r
length(table(tx2gene$TXNAME))
length(table(tx2gene$GENEID))
```

Honestly there is something weird (possibly the new genome annotation) as there should only be about 15-20k genes and 35K transcripts.. So....

Now we can go ahead and read in the input — this will automatically sum results to the gene level.  You can check out the tximport documentation for some other, potentially useful options.


```r
txi <- tximport(quant_files,
  type = "salmon",
  tx2gene = tx2gene)

# Let's also take a quick look at txi

summary(txi) # why 274560?
str(txi)
head(txi$counts) # note these are not integers!


# We can look at patterns of correlations among samples
cor(txi$counts)


pairs(txi$counts[,1:2], 
      pch = 20)
```

## Setting up our experimental design.

DESeq2 needs the information about your experiment set up, so it knows the various predictors in the model (in this case genotype and background). The easiest way to do this is by setting it up as a data frame in R (which is a specialized version of a list). I will (time permitting) show you a more general way of doing this with another example, but for now we are explictly writing this out.


```r
tissue <- c(rep("genital", 5),
            rep("wing", 6),
            rep("genital", 3),
            rep("wing", 6))

tissue <- as.factor(tissue)
length(tissue)

temperature <- c(rep(17, 11), rep(24,9))
temperature <- as.factor(temperature)
length(temperature)

food <- c(rep("fed", 3),
          rep("starved", 2),
          rep("fed", 3),
          rep("starved", 3),
          rep("fed", 2),
          rep("starved", 1),
          rep("fed", 3),
          rep("starved", 3))
food <- as.factor(food)
length(food)

lane <- c(2,4,5,5,2,3,2,4,4,3,2,4,2,3,2,4,3,3,2,4)
lane <- factor(lane) # we will want to treat this as a factor
length(lane)

rna.design <- data.frame(sample=samples,
  file=quant_files,
  tissue=tissue,
  food=food,
  temperature=temperature,
  lane = lane)

rna.design

# and we can start with a simple model (back to this later)
load.model <- formula(~ tissue)
```

Below is the crucial function, `DESeqDataSetFromTximport()`, that gives you a DESeq2 count matrix from the txt object. interestingly, this is the step that converts the estimated counts to integers (again we can take a look at the quant_files).


```r
all.data <- DESeqDataSetFromTximport(txi, 
                                     rna.design,
                                     design=load.model)
```

## Data is in, now what?

This is normally a good opportunity to do some simple visulizations to look at the distributions of the estimates and the correlations among samples (what should we be looking for). 


## Preliminary Quality Control analysis
Before we begin any real analysis. It pays to take some looks at the data. I am not going to go through a full exploratory data analysis session here. But some obvious plots

It is well known that there can be substantial lane to lane variation. For this experiment, it was designed so that a number of samples were run in each lane (barcoded), in a randomized design. This enables us to control for lane effects if necessary. As it turns out I picked a somewhat useless sub-sample of the full data set, so we can not look at the lane effects (as we don't have enough samples in each lane for this subset of data we provide). But normally do something like this (and include a lane effect at a covariate)

First we create a DESeq data object using our counts, experimental design and a simple statistical model (more on this later)


```r
load.model <- formula(~ lane)

test_lane_effects <- DESeqDataSetFromTximport(txi,
  rna.design, design=load.model)

test_lane_effects2 <- DESeq(test_lane_effects)
# We now fit the simple model
```

This generates a fairly complex object

```r
str(test_lane_effects2)
```

For the moment we can ask whether any genes show evidence of different expression based solely on lane to lane variation. We use the results() to summarize some of the results.


```r
test_lane_effects2_results <- results(test_lane_effects2, alpha = 0.05)
# alpha = 0.05 is the  "cut-off" for significance (not really - I will discuss).

summary(test_lane_effects2_results)
# 2 genes which may show  evidence of lane effects, but this is a bit incomplete for the full data set.

head(test_lane_effects2_results)

# let's re-order the data to look at the two genes.
test_lane_effects2_results <- test_lane_effects2_results[order(test_lane_effects2_results$padj),]
```

We can also plot the mean-dispersion relationship for this data.


```r
plotDispEsts(test_lane_effects2)
```
Let's talk about what this means.

### Principal Components analysis and hierarchical clustering are useful tools to visualize patterns (and to identify potential confounds)

 We can also use some multivariate approaches to look at variation. For PCA (checking it with a "blind" dispersion estimate to look for any funky effects. Not for biological inference).


```r
for_pca <- rlog(test_lane_effects2, 
                blind = TRUE)
dim(for_pca)
```
`rlog` is one approach to adjusting for both library size and dispersion among samples. `blind=TRUE`, has it ignore information from the model (in this case lane).


```r
plotPCA(for_pca, 
        intgroup=c("lane"),
        ntop = 5000) 
```

The `plotPCA()` function is actually just a wrapper for one of the built in functions for performing a principle components analysis. The goal of this (without getting into the details for the moment) is to find statistically independent (orthogonal) axes of overall variation. PC1 accounts for the greatest amount of overall variation among samples, PC2 is statistically independent of PC1 and accounts for the second largest amount of variation. By default the `plotPCA` function only plots the first couple of Principle components. In this case it explains just under 80% of all of the variation among the samples. However, I highly recommend looking at the plots for higher PCs as well, as sometimes there is something going on, even if it only accounts for a few % of variation.

If you want to see what this wrapper is doing we can ask about this particular function


```r
getMethod("plotPCA","DESeqTransform")

# or 
DESeq2:::plotPCA.DESeqTransform
```
It is very easy to modify this function if you need to.

By default this only examine the top 500 genes. Let's look at 2000 or more to get a sense of the stability of the pattern.


```r
plotPCA(for_pca, ntop = 1000, intgroup=c("lane")) 
```


### Back to the analysis
While there is some lane effects based on our initial statistical check, and visual inspection of PC1 VS. PC2. However, there clearly is a pattern here, and it has nothing to do with lane. 

We can quickly take a look to see if this pattern shows anything interesting for our biological interest. However, this is putting the cart before the horse, so be warned.


```r
plotPCA(for_pca, ntop = 5000,
        intgroup=c("tissue", "food", "temperature"))
```

Not entirely clear patterns of clustering here. Play with this changing the number of genes used. Also keep in mind that it is by default only showing the first two principal components of many, so it may not be giving a very clear picture!


### We can also use some hierarchical clustering to further check for lane effects or for clustering based

For distance matrix for clustering QC

```r
rlogMat <- assay(for_pca) # just making a matrix of the counts that have been corrected for over-dispersion in a "blind" fashion
distsRL <- dist(t(rlogMat)) # Computes a distance matrix (Euclidian Distance)
mat <- as.matrix(distsRL)  # Make sure it is a matrix
```

We need to rename our new matrix of distances based on the samples.

```r
rownames(mat) <- colnames(mat) <-   with(colData(test_lane_effects2), paste(genotype, background, sep=" : "))

hc <- hclust(distsRL)  # performs hierarchical clustering
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)  # picking our colours
```

Now we generate the plot

```r
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm = TRUE, trace="none",
          col = rev(hmcol), margin=c(8,8))
```

 While we checked that there was modest evidence of lane effects , we could keep it in the model as it seems to have little effects. However, it is using up "df" (in the parameter estimation sense), so it may be worth ultimately getting rid of it.



## Proceeding with the real analysis we care about!
Given the results from above, I am removing lane entirely.

Let's fit and run the model.

```r
load.model <- formula(~ lane + tissue)

test_tissue_effects <- DESeqDataSetFromTximport(txi,
  rna.design, design=load.model)

test_tissue_effects2 <- DESeq(test_tissue_effects)
```

Now we can look at some of the results. First we recheck the dispersion estimates estimated in our model


```r
plotDispEsts(test_tissue_effects2)
```
Not much different. A few outliers though, and it may be worth later going in and checking these.

Let's get going with things we are interested in, like looking for differentially expressed genes across genotypes. We can start by doing a visualization using a so-called MA plot (look it up on wikipedia, then will talk)


```r
plotMA(test_tissue_effects2, ylim =c(-6, 6))
```
 A few things to note. The points coloured in red are the genes that show evidence of differential expression. The triangles are ones whose log2 fold change is greater than or less than 1 ( i.e. 2 fold difference). Please keep in mind that many genes that have small fold changes can still still be differentially expressed. Don't get overly hung up on either just fold changes alone or a threshold for significance. Look carefully at results to understand them. This takes lots of time!!!
 
Let's actually look at our results. DESeq2 has some functions to make this a bit easier, as the object we generated above (of class `DESeqDataSet`) is quite complex.


```r
tissue_results <- results(test_tissue_effects2, 
                            alpha = 0.05)
print(tissue_results)

head(tissue_results)

summary(tissue_results)
```

A few things to note. By default alpha is set at 0.1 This is a pretty liberal "threshold" for assessing "significance". While it is a much larger conversation, I do not recommend getting hung up too much on statistical significance (since you have estimates of effects and estimates of uncertainty). In particular if you have specific questions about specific sets of genes, and if these were planned comparisons (i.e. before you started to analyze your data) you should focus on those.

You will see that DESeq2 also will throw out genes that it deems outliers or have very low counts. It does this so that when it is correcting for multiple comparisons it does not need to include genes that should not be analyzed (as the counts are too few.)


Let's take a look at the results.


```r
# reorder
tissue_results <- tissue_results[order(tissue_results$padj),]

tissue_results[1:10,]
```

What do you note about the differences. Please compare these different lists of genes.

## More complex models

While everything is stored, by default DESeq2 is printing its evaluation of the final term in the model. We can look at these in the model we actually want to fit (unlike the naively simple model above)

Let's start by examining the effects of genotype (like we did above), but by first taking into account difference in different wild type genetic background and any residual lane effects



```r
load.model <- formula(~ lane + tissue + food) # Let me go over the matrix rank issue here

test_G <- DESeqDataSetFromTximport(txi,
  rna.design, design=load.model)

test_G_2 <- DESeq(test_G)

plotDispEsts(test_G_2)

plotMA(test_G_2, ylim =c(-4, 4))

G_results <- results(test_G_2, alpha = 0.05, pAdjustMethod="BH")
summary(G_results)

# reorder
G_results <- G_results[order(G_results$padj),]
G_results[1:10,]
```

### Interaction terms in models. 
In our case we are as much interested in genes that show an interaction between genotype and background as those that show just an effect of genotype. However for the subset of all of the samples we are looking at, we are going to run into estimation issues.


```r
load.model <- formula(~ lane + tissue + food + tissue:food) # Let me go over the matrix rank issue here

test_FxT_effects <- DESeqDataSetFromTximport(txi,
  rna.design, design=load.model)

test_FxT_effects2 <- DESeq(test_FxT_effects)

plotDispEsts(test_FxT_effects2)

plotMA(test_FxT_effects2, ylim =c(-4, 4))

FxT_results <- results(test_FxT_effects2, alpha = 0.05, pAdjustMethod="BH")
summary(FxT_results)

# reorder
FxT_results <- FxT_results[order(FxT_results$padj),]
FxT_results[1:14,]
```

Relatively few "significant" genes. Just to keep in mind. We expect a priori, with no true "significant" hits an approximately uniform distribution that looks like (for 13140 genes)


```r
p_fake <- rbeta(13140, 1,1) # you could also use runif(12627,1,1)
hist(p_fake)
```

But we actually observe

```r
hist(FxT_results$pvalue)
```

False Discovery Methods (FDR) exploit this.
