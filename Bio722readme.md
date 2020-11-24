# BIO722 Fundamentals of genomic data analysis

Winter/Spring 2021 (expected to be taught love online via zoom or teams)

## Instructors

Dr. Ben Evans

Dr. Brian Golding

Dr. Ian Dworkin (dworkin "at" mcmaster "dot" ca)

## Course outline

This course will help graduate students further develop their computational skills with respect to genomic analysis. The goal is very much on the practical side of the computational data analyses for genomic data (RNAseq, DNA resequencing experiments, variant calling, gene ontology analysis among others). The first third to one half of the semester will be led by the instructors. Classes will be a mix of lecture style, discussion and tutorial activities. The second half of the semester will be student led workshops style lecture/tutorials on topics chosen by students (and approved by instructors).


## Course goals

1. To provide a practical (approximately 2/3) and conceptual introduction (1/3) to fundamental skills and background for genomic analysis. This includes QC, analysis, and visualization steps among others.

2. To facilitate students to independently develop their computational skills and background for genomic analysis.

3. Help students to develop appropriate habits for complete workflows, with appropriate data "cleaning" and QC steps all along the pipeline. How to build in sanity checks into your pipeline.

4. Help students develop skills and need for reproducible research. We will be using git (and github) for version control. Students will be expected to develop and maintain their own git repositories for this class that will be shared with other class members and instructors.

## Course topics

1. RNAseq: The first 3-4 weeks will be led by Dr. Dworkin, and will go over conceptual and practical considerations for RNAseq experiments for differential expression. (experimental design, QC, mapping, read counting, more QC, statistical modeling, regularization, dealing with multiple comparisons, did I mention QC...). We will use a combination of mapping software (Salmon, STAR) and R libraries for differential expression (DEseq2, limma/voom, edgeR) along with a number of tools to perform quality control throughout the pipelines. Visualizations will be done through plotting using R.

2. Burrows Wheeler Transform (Dr. Golding). A conceptual overview of the Burrows Wheeler transform which is fundamental for many of the mapping tools that are commonly used.

3. DNA "resequencing" experiments, syntenic mapping, reduced representation libraries and variant calling. Dr. Evans will lead 1-2 lecture/tutorials on using DNA resequencing experiments with the goal of mapping and identifying variants (SNPs, indels etc) and using these for downstream analyses (like population genetics, possibly using PLINK for genome wide association mapping). Dealing with range data might be introduced here (or possibly earlier)

4. The remaining subjects will be jointly determined by the students and the instructors and will be presented/facilitated by student groups. Topics may include things like gene ontology analysis, transcriptome assembly, genome assembly, CHiPseq, ATACseq (or related approaches like FAIREseq, DNAseIseq etc), identifying novel transcripts (and transcript level abundance analysis), genome annotation...

## What we do not cover in this course

- We do not review the fundamentals of programming in UNIX, R, python, etc. We assume students have previous background or taken a course like Bio720 or equivalent.

- We do not discuss or do any wet lab experiments. So we do not go into detail on DNA, RNA, protein (etc) extraction, library prep etc.

- We do not explicitly go into an overview of all types of genomic analyses for different data types. So some applications may not be discussed (but could still potentially be used for independent projects).

- We do not go into any details of all sequencing technologies. Most examples provided by Instructors will likely be from Illumina short read sequencing experiments. However, if students so choose, group presentations/tutorials can use long read examples  (or other approaches).

- While we might use some simulations approaches (i.e. simulating whole genomes or transcriptomes) the main focus of the course is on the data analysis of genomic data *per se*.

## Assumed background for students taking this class

- Unlike [Bio720](./README.md) this class assumes you already have some familiarity with both UNIX shell scripting and at least one other scripting language (R, python, perl, etc). Activities in the class will provide some example code (mostly UNIX and R), but students will be expected to be writing much of their own code and examples.

- If you have not taken Bio720, but have a reasonable background in both UNIX shell scripting and a standard scripting language, please contact the Instructors.


## Course evaluation

The grade for this course will be based on a Final independent Project and on student led/facilitated workshop style tutorials (short lecture and practical). The exact number and format of student led tutorials will be determined the first week of class (once we know enrolment).

 While there will be practical activities during the first few weeks of the class , these will not be formally marked. However, if you want the instructors to check over your code.
