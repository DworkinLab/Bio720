# Bio720


## Course Summary
This is my (ID) github repo for students for  Bio720 (introduction to computational skills for biologists).

[Link to R portal for class](https://github.com/DworkinLab/Bio720/blob/master/Introduction_to_R.md)

[Link to UNIX and shell portal for class](http://helix.mcmaster.ca/720.html). Will be updated for 2nd half of class.

## Instructors
[Dr. Brian Golding](http://helix.biology.mcmaster.ca/)


[Dr. Ian Dworkin](https://dworkinlab.github.io/)

## Class times and locations

BIO720, "lecture" Mondays, 9:30-10:20, ETB 227 **PLEASE NOTE ROOM CHANGE, starting September 16th**  

BIO720, "practical lab" Tuesdays, 9:30 - 11:20 ETB228

## Information on submitting assignments and grading

[Please follow this link](./BIO720_assignments_HowToSubmitGrade.md)


## Background assumed for students for Bio720
For this class we are not assuming students have background in programming/scripting, nor in bioinformatics. We do assume that students have a working knowledge of basic molecular biology and genetics and have basic familiarity with using computers. i.e. you can figure out how to install basic software on a Mac (OS X) or a PC (Windows).

## Background assumed for students for BIO722 (when offered)

That you have taken BIO720 or equivalent. Are comfortable at the UNIX shell, with at least one scripting language and able to learn a second. Know how to remotely access computers (ssh etc).

# What students will need
A laptop with internet access and the ability to install several programs (in particular, `R`, `Python` and a shell emulator if not using a Mac (OS X) or linux. If not using a Mac, we highly recommend you look into making your computer a dual boot machine (so it can boot both in Windows and UNIX)

## Course goals
The primary goal of this course it to provide graduate students an opportunity to develop fundamental computational skills necessary to go on and (in the future) develop the appropriate (and more advanced) skills for bioinformatics, genomics, etc.

## **What this course is not**
Because of limitations of time, we are purposefully making this a course about fundamental skills. As such, this course will not cover in any detail:

    - Genomic analysis pipelines (RNAseq, variant calling and populations genomics). These are covered in BIO722.

    - Theory of computer science (nor theory on programming, algorithms, data structures etc).

    - Despite using `R` for much of this course, it is most definitely not a statistics course. Bio708 (taught by Dr. Ben Bolker and Dr. Jonathan Dushoff and Dr.Dworkin) is such a course (also using R as the primary programming environmental for statistical modeling.)

    - A bioinformatics class (i.e. we will not teach any conceptual or theoretical background in bioinformatics. All examples will be real examples, but mostly to illustrate the computational skills necessary to run an analysis, not the why).


## Learning Objectives

### Topics (some TBD)
 It is important to note that in order to keep things flexible depending on how things go with the class, these topics are subject to change if necessary. We will discuss in class.

A. Fundamentals of programming using `R` (Ian). [Link to R portal for class](https://github.com/DworkinLab/Bio720/blob/master/Introduction_to_R.md)
  1. Fundamentals of programming in R.
  2. How to avoid repetitive strain injury while programming. Control flow in R (`for` loops, `if else`, etc). Using the `apply` family of functions in `R`. Simple  simulations.
  3. Working with data in `R`. Getting data in. Data munging (subsetting, merging, cleaning). Working with strings in `R`.
  4. Basics of plotting in `R`. Other topics TBD.
  5. Reproducible research using markdown for reports and git for version control.
  6. An Introduction to bioinformatic tools in R. Primarily an introduction to BioConductor, and genomic range data.

B. Introduction to UNIX and the command line. (Brian)
  1. Introduction to basic shell commands, logging onto remote systems
  2. Standard UNIX utilities that make your day to day computer work (and bioinformatics) easier.
  3. using pipes in UNIX (and the model of streaming data), batch processing of data.
  4. Writing shell tools.
  5. Using your UNIX skills for practical bioinformatic problems (probably setting up a BLAST database, and querying some sequences)
  6. (maybe) Regular expressions are you friend. No really. Using `grep` and its variants (i.e. `agrep`) and `sed` and `awk` for file manipulation and processing

### Learning outcomes

After successfully completing this course you will:

    - Have a much higher degree of comfort using your computer!

    - Be able to write custom UNIX shell scripts to do file copying, moving, editing, parsing and manipulation.

    - Be able to write simple R programs to do simple simulations, data parsing (munging), plotting.

    - Be able to perform computationally reproducible research, and use version control on your source code.

    - Be able to utilize genomic range data and incorporate simple genomic features.

    - Understand the fundamental framework of UNIX programs, scripting and why streaming data is so useful for genomics and bioinformatics.

    - Know that troubleshooting for installing and using programs, and troubleshooting when writing and using code are normal. You will have developed some tenacity in dealing with such issues and have some ideas on how to approach finding solutions (including your *google-fu*).

## Recommended books.
You are responsible for ordering your own copies of these books if you want physical copies. McMaster has access to both books online. Both are excellent with some  overlap. 

[Bioinformatics Data Skills, **BDS** by Vince Buffalo](https://mcmaster.primo.exlibrisgroup.com/permalink/01OCUL_MU/deno1h/alma991032459102207371). **HIGHLY RECOMMENDED, and we have free access through the McMaster library** This book fills an important gap in that is oriented towards the day to day skills for anyone working in the fields of genomics and bioinformatics. In addition to covering the basic UNIX skills (and why we use UNIX in bioinformatics and genomics), it also covers subjects like overviews of the essential file types (`.fasta`, `.fastq`, `.gff`, etc) that are ubiquitous in the field. There is also a nice, but brief introduction to the essentials of `R`, using bioconductor and in particular range data, and two important chapters on how to organizing (and maximize reproducibility) of computational projects. The author wrote this in their first year of graduate school, so definitely worth supporting them!

[Computing Skills for Biologists: A toolbox, by Stefano Allesina and Madlen Wilmes](https://mcmaster.primo.exlibrisgroup.com/permalink/01OCUL_MU/deno1h/alma991032955925007371). **HIGHLY RECOMMENDED, and we have free access through the McMaster library**. This book provides a nice, gentle introduction to the basic computational skills all biologists should have. In particular, with introduction to using the UNIX command line, shell scripting, basic `python` programming, regular expressions, working on remote machines and a few other topics. The book is written to be agnostic with respect to discipline (i.e. it is not a bioinformatics book *per se*), but does a great job of being both very accessible and immediately useful. If you plan to continue in computational research, this is a fantastic resource.


## Important websites

[R tutorials and screencasts](https://github.com/DworkinLab/Bio720/blob/master/Introduction_to_R.md). A link to the exercises, in-class activities, playlists for screencasts I have put together for the `R` tutorials. I will also be putting assignments up here. I will be adding more as the semester progresses. Mostly we will be using the excellent [Datacamp](https://www.datacamp.com/) online interactive 'courses' for the introductory stuff, and moving on from there.

[For Brian's section](http://helix.mcmaster.ca/720.html). This will have pertinent links to Brian's section of the course.

## Readings
Week 1 - BDS Chapter 1, Chapter 2 pages 21-30, Chapter 3 pages 37-45, Chapter 4 pages 57-59.

Also check out [here](https://github.com/DworkinLab/Bio720/blob/master/IntroductionMarkdownAndVersionControl/Bio720_IntroductionMarkdown.md#a-few-words-on-project-organization) for a review of organizing computational data analysis projects. You don't need to read about version control or using markdown **yet**.

Week 2 (Beginning R) - Start with DataCamp assignments (Courses **Introduction To R** and **Intermediate R**). Associated readings in BDS chapter 8 pages 175-206.

Week 3 - Chapter 3 pages 45 - 56, Chapter 7 125 - 156. Maybe also worth looking at pages 395 - 398 in Chapter 12.

Week 4 - Chapter 7 140-145, 157-169 might be really useful. I also recommend the first tutorial on regular expressions listed [here](https://github.com/DworkinLab/Bio720/blob/master/Some_RegEx_resources.md). This takes you gently through regular expressions and within an hour you will realize what amazing things you can do.

Week 5 - For more information on some of the basic file types used in genomics (.fastq, .SAM, .BAM) see chapter 10 (only 13 pages) and chapter 11 (pages 355-365) in BDS. I also suggest reading chapter 6 (pages 109-115).


## Weekly materials (in class)

Week 1:

Welcome to BIO720 - [slides](./lectures/) or [markdown](./lectures/Welcome_BIO720.md)

Short introduction to reproducible research - [slides](./lectures/) or [markdown](./lectures/BIO720_VersionControl_reproducibleResearch.md)

Time permitting some activities introducing R.

Week 2:
**Monday**

- Pseudo-code.

- Working with data in R using your new skills.

  -- [You will want to be able to access this script to download the RNAseq data we will be playing with](./R_exercises/R_activities_InClass_Sept9_2024.Rmd)
  
Week 3:  

**PLEASE NOTE THE ROOM CHANGE for Monday, ETB 227** 

- Monday and Tuesday (control flow)
    - [control flow fundamentals](./R_exercises/Bio720_R_InClass_Control_flow_fundamentals.Rmd)
    
    - [This is the script we will be using in class](./R_exercises/Bio720_R_InClass_Control_flow_worked_example.Rmd)
        - [The script with answers](./R_exercises/Bio720_R_InClass_Control_flow_worked_answers.Rmd)

Week 4: 

- Monday (making and breaking functions, in class activity) 
    - [The script to start with](./R_exercises/MakingAndBreakingFunctions_InClass.Rmd)
    
- Tuesday, working with strings in `R`.
    - [A short introduction to working with strings](./R_exercises/Bio720_AMorningofStrings_short_InClass.Rmd)
    
        - [answers: A short introduction to working with strings](./R_exercises/Bio720_AMorningofStrings_short.Rmd)
    
    - [Data munging in base R](./R_exercises/Bio720_R_week3_DataMunging_InClass.Rmd)


## Assignments

[Assignment 1, due Monday, September 23rd, 2024](./assignments/Bio720_2024_R_Assignment_1.pdf)
