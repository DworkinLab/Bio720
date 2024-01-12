---
title: "BIO722 - Genomic Analysis"
author: "Instructors: Ian Dworkin, Ben Evans, Brian Golding"
date: "12 Jan 2024"
output:
  ioslides_presentation: 
    fig_height: 4
    fig_retina: 1
    fig_width: 6
    keep_md: yes
  slidy_presentation: 
    fig_retina: 1
    fig_width: 6
    keep_md: yes
  beamer_presentation:
    incremental: no
editor_options:
  chunk_output_type: console
---



## Outline for today
- Introductions around the class
- What to expect in this class
- Course outline
- Course projects and evaluation
- Discussion: Why do we use genomics? What is good for? What are its limitations?

## What to expect in this class
- This class is very much what you make of it.
- While we will be running lectures and tutorials for the first half of the class, it is mostly independent and small group project work.
- Genomic analysis needs to be learned through actively doing it, discussing the work, and refining your thinking and analyses.

## Semester overview
- Introduction to differential gene expression analysis using RNAseq (and fundamentals of genomic analysis). ~3 weeks.
 - Syntenic mapping of sequence reads to transcriptomes and genomes
 - Counting and analysis of mapped reads
 - Differential Gene Expression analysis of RNAseq data
 - Quality assessment at multiple stages of analysis

- Project pitches by students (first week of February)

- Brian (TBD) - Likely some advanced UNIX tools to make large scale analyses run smoothly (2 weeks)

- Ben Evans - Population Genomics, variant calling, GWAS (TBD)

- The rest of the semester will be group workshop/tutorials presented to the class. 

## Course evaluation breakdown
Project Proposal: 10%
 - First week of February
 - Written (2 pages maximum)
 - Class pitch (10 minutes + 5 minutes questions)
 
Group Project: 30%
 - workshop/tutorial led by small groups (2-3 people)
 - Throughout March

Progress report for independent project: 5%
 - Early March
 - Written, 1-2 pages
 - We may do a short in class oral update

Final Project: 45%
 - Due mid-late April
 - Written paper and annotated and reproducible code (github)

Class Participation: 10%

## Possible topics for group workshops/tutorials
- Transcriptome Assembly
- Differential transcript (or exon level) analyses (DTE, DTU, etc)
- meta-genomics
- DNA methylation (analysis of bisulfite sequencing)
- Chromatin accessibility (ATACseq, DNAse1seq, FAIRE, etc)
- Differential co-expression analysis
- How to combine multiple genomic analyses together (multi-omics)
- Genome annotation
- Gene Ontology analyses


## Why genomics
- Before we spend the rest of the semester just "doing it" with respect to genomic analysis, let's spend a bit of time on
 - **Why** do we do genomics?
 - **What** can (in an optimal world), genomic studies tell us?
 - What is it (in a real world) likely to be able to tell us?
 - What can't genomics tell us?
 - Are genomic analyses hypothesis generating? Testing? Both? Neither?
 - **How** to (not) lie with genomics.
