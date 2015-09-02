
The current plan is ...

Four lectures each:

BG: The basic skill set for practical bioinformatics research (with NGS data)

1. UNIX cmd line (2 classes): redirection, pipelines, etc, Editor - nano
2. Programming in R *(maybe to ID?)*
3. File formats - fastq, gff, (Quality exploration - fastqc), Trimming reads (Trimmomatic(
4. Hash tables, SA arrays, BWT, FM datastructures
   syntenic aligners - sam file format (bwa or bowtie2)


BE: Identifying and utilizing polymorphism data from Illumina resequencing (and reduced representation sequencing).

Introduction to the Macaque data
Data types for identifying polymorphism (reduced representation approaches in particular).
SNP calling/genotyping: GATK (and base recalibration, deduplication, unified genotyper VS. haplotype caller), vcftools, vcf file formats, bed file format

ID: Introduction to using RNAseq for the analysis of differentially expressed genes

0. Markdown (this probably belongs with R) 
1. Transcriptome Assembly using Trinity (ID can do the practical part: Maybe BG for theory?), and how to interpret and use the **de novo* transcriptome assembly. 
2. Considerations for types of features to use (genes, transcripts, kmers), How to count for RNAseq (what features to use, different approaches to deal with multi-mapped reads)
3. Basic Statistics for the analysis of count data for differential expression
4. Differential expression analysis. Examining data for technical issues. Dealing with multiple comparisons (FDR). 
    
