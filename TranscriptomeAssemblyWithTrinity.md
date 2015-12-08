# Transcriptome assembly with Trinity.

For this tutorial we will just go through the code, and not run it (it takes too long for our time).

In addition to this tutorial, I recommend a number of other resources.

First, at the NGS summer workshop, we have done transcriptome assembly a number of times, with different instructors, data sets and tools. This will give  you a tiny flavour for the diversity of approaches.

-  [Matt MacManes transcriptome evaluation tutorial (2015)](http://angus.readthedocs.org/en/2015/MacManes_Trinity.html)

- [Meg Staton's tutorial on using transrate](https://github.com/ngs-docs/angus/blob/2015/transrate.rst)

- [Titus Brown's "complete" protocol](http://angus.readthedocs.org/en/2014/eel-pond.html). This includes using his (CTBs) diginorm tool as opposed to the digital normalization in Trinity.

- [An older tutorial for velvet-oases] (http://angus.readthedocs.org/en/2013/transcriptome_de_novo_assembly.html)

- [Main Trinity tutorial](https://github.com/trinityrnaseq/RNASeq_Trinity_Tuxedo_Workshop/wiki). There is loads of info on the trinity page, and they do their best to guide you step-by-step. Also see the protocol in [Nature Protocols](http://www.nature.com/nprot/journal/v8/n8/pdf/nprot.2013.084.pdf)

## What steps will we do.
1. Quality & Adapter trimming using Trimmomatic
2.

## 1. Quality Trimming with Trimmomatic
When we did quality trimming with Trimmomatic for the RNA-seq reads for *counting* we mostly made sure to trim off the adapter sequences, and VERY low quality reads. However, as sequencing error (like polymorphism) may introduce issues during assembly resulting in false contigs, I generally recommend somewhat more severe trimming for *de novo* transcriptome assembly. It is worth noting that [not everyone recommends this](https://impactstory.org/MatthewMacmanes/product/7v4g4tw0xryyesqi5m4o6og7/fulltext). Instead Matt MacManes suggests [error correction](https://peerj.com/articles/113/).

```bash
java -jar /usr/local/trimmomatic/trimmomatic-0.33.jar PE -threads 8 \
${dir}${name}_R1_001.fastq.gz ${dir}${name}_R2_001.fastq.gz ${dir_out}${name}_1_pe.fastq\ ${dir_out}${name}_1_se.fastq ${dir_out}${name}_2_pe.fastq ${dir_out}${name}_2_se.fastq\ ILLUMINACLIP:/usr/local/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10\
LEADING:3 TRAILING:3 MAXINFO:20:0.6 MINLEN:36
```

So the big difference is the maxinfo setting (0.6) which aims for moderate-high degrees of trimming (as opposed to 0.1-0.2 which is more liberal about quality). The `MINLEN` means that any read that is cut to less than 36bp will be removed entirely. As always you would write a `for` loop to run through this for pairs of samples, or use [parallel](https://www.gnu.org/software/parallel/parallel_tutorial.html) or the clusters job scheduler to run them.

At the end of this you will (for each sample) end up with a pair of files for the remaining paired end sequences, and then two files with the remaining single end sequences (where one of the two reads did not survive quality filtering). Whether you want to keep the single end remaining reads (across samples) is upto you. Certainly, `Trinity` can handle them. I recommend using `cat` to concatenate all of the single end (SE) fastq files.

i.e.:

```bash
cat *_se.fastq > all_se.fastq
```

This will make your life a bit simpler. If you do this, it is recommended that these single end reads are added to the last paired end (PE) file of the samples going into trinity (again using `cat`). So this may look something like this:

```bash
cat name_of_last_pe.fastq all_se.fastq > name_of_last_pe.fastq
```

Note I **did not** change the name of the last file. This is because it is one of two files (paired end) that need to match in name. See the Trinity manual for more details.


## 2. Running Trinity

Since Trinity can handle the digital normalization, it is fairly easy to use. As I mentioned in lecture, Trinity seems to have fairly consistent results to changes in parameter values, and even the vanilla settings seem to work reasonably well for many researchers.

```bash
/usr/local/trinity/Trinity --seqType fq --normalize_reads --JM 20G --CPU 8\
 --left M160_lg_male_genitalia_AGTTCC_L004_1_pe.fastq,M160_lg_male_hdhorn_TTAGGC_L005_1_pe.fastq,M160_lg_male_thxhorn_GTGAAA_L002_1_pe.fastq,M160_lg_male_wings_ATGTCA_L007_1_pe_se.fastq\
  --right M160_lg_male_genitalia_AGTTCC_L004_2_pe.fastq,M160_lg_male_hdhorn_TTAGGC_L005_2_pe.fastq,M160_lg_male_thxhorn_GTGAAA_L002_2_pe.fastq,M160_lg_male_wings_ATGTCA_L007_2_pe.fastq\
   --min_contig_length 300 --output ./trinity_out
```

Let's unpack this a bit
- `seqType` is just the sequencing file type (in this case fastq)
- `normalize _reads` Whether to use Trinity's built in digital normalization.
- `JM` Settings for how much memory can be used. Multiply this by number of processors to get RAM used (in this case 160GB of RAM)
- `CPU` how many threads you can use
- `left` all of the files associated with the "left" hand reads of the pairs. Ditto for right.
- `min_contig_length` the minimum accepted contig length for output (shortest transcripts you will expect.)

Since this can take a LONG TIME to run. For a single individual beetle (across tissues) it was typically 18-20 hours. For the multi-individual runs, it was 3-5 days. Often it will crash when running. Thankfully Trinity is pretty good about saving intermediate files, and if you look at the help, they can show you how to run it from the last checkpoint.

Another argument to keep in mind is `--jaccard_clip` which should reduce the number of fused (chimeric) transcripts, but at a cost of speed (often doubling or tripling time needed). However, they suggest it can reduce this by ~50%, so it may well be worth it.


## 3. Post-assembly
So you know have an assembly. IN this case you may get output (`Trinity.fasta`) with all of your contigs (both components, which represent "genes" and transcripts within the components). The authors of Trinity have made a number of very helpful tools to evaluate this. To get some basic summary stats you can use;

```bash
/usr/local/trinity/util/TrinityStats.pl ~/yourDirectory/Trinity.fasta
```


Which should give you output that looks like:
```bash
################################
## Counts of transcripts, etc.
################################
Total trinity 'genes': 20014
Total trinity transcripts: 34802
Percent GC: 36.49

########################################
Stats based on ALL transcript contigs:
########################################

Contig N10: 6855
Contig N20: 5024
Contig N30: 4002
Contig N40: 3272
Contig N50: 2676

Median contig length: 1118
Average contig: 1687.61
Total assembled bases: 58732247

#####################################################
## Stats based on ONLY LONGEST ISOFORM per 'GENE':
#####################################################

Contig N10: 6638
Contig N20: 4879
Contig N30: 3847
Contig N40: 3079
Contig N50: 2507

Median contig length: 869
Average contig: 1485.75
Total assembled bases: 29735717
```

They also have nice tutorials to help you evaluate some basic quality aspects of your transcriptome, using blast. See [here](http://trinityrnaseq.github.io/analysis/full_length_transcript_analysis.html) for more details.

Generally this involves making a protein database:

```bash
/usr/local/blast/2.2.29/bin/makeblastdb -in uniprot_sprot.fasta -dbtype prot
```

Where I downloaded the most recent copy of uniprot and made a blastable database for it.

I then used blast against the transcriptome I just assembled

```bash
/usr/local/blast/2.2.29/bin/blastx -query ../yourDirectory/Trinity.fasta -db uniprot_sprot.fasta -out blastx.outfmt6 -evalue 1e-20 -num_threads 6 -max_target_seqs 1 -outfmt 6
```

In this case I am making an output `blastx.outfmt6` where I am specifying only the tophit (`max_target_seqs`) with a minimum evalue and a specific output format (`outfmt`).

I then used one of the helper Trinity perl scripts to assess how well my de novo assembled contigs blasted against the database

```bash
/usr/local/trinity/util/analyze_blastPlus_topHit_coverage.pl ~/yourDirectory/blastx.outfmt6 ~/yourDirectory/Trinity.fasta ~/databases/uniprot_sprot.fasta
```

(note you may not have a single folder for your blastable databases, so you may need to change the path).

This should give you some output like this:

```bash
100 2570 2570
90 944 3514
80 649 4163
70 509 4672
60 472 5144
50 525 5669
40 484 6153
30 477 6630
20 391 7021
10 111 7132
```

So we have 2570 transcripts that have 100% identity (over some pre-specified length), 3514 with at least 90% identity and so on. You can then do some scripting to investigate what they are. Some folks also like using blast2go for these steps as well.

## 4. What's next

If you are like me, you will spend a fair bit of time trying to assess the quality of the transcriptome (using something like transrate), and also trying to simplify the transcriptome using CAP3 (contig assembly program) or CD-hit-est. In my experience how well either of these do depends a great deal on the parameters you set, and I do not feel like enough of an expert to make firm suggestions (but I think Brian has used CD-hit-est a fair bit).


I have recently been working with [TransPS](https://bioinformatics.cs.vt.edu/zhanglab/transps/) which aligns the contigs to a genome or transcriptome (a high quality one), queries are then ordered based on the quality of the alignment and non-redundant contigs are merged into a contig. I have not yet done enough evaluation to assess how well this works relative to other approaches yet. Currently I only have TransPS installed locally on Brian's cluster and it took a fair bit of futzing to get it to compile. However, once you do  running it is pretty straightforward. Note I re-ran blast and allowed for more matching hits.

```bash

/usr/local/blast/2.2.29/bin/blastx -query rinity.fasta -db ~/databases/Tc_NCBI_protein.fa -out Tc_blastx_20Hits.outfmt6 -evalue 1e-20 -num_threads 16 -max_target_seqs 20 -outfmt 6

perl ~/temp/TransPS1.1.0/transps.pl -t combined_trinity_out.Trinity.fasta -b Tc_blastx_20Hits.outfmt6
```
