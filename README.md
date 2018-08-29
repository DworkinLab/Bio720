# Bio720

## Course Summary
This is the page for Bio720 (2016), a practical introduction to fundamental computational skills for biologists. This is taught through the Biology Department at McMaster University, but most course materials are freely available to anyone interested.

## Instructors
[Dr. Brian Golding](http://helix.biology.mcmaster.ca/)

[Dr. Ian Dworkin](https://scholar.google.com/citations?user=Iium3AEAAAAJ&hl=)

## Class time and location
TBD. First organizational meeting will be held in the Life Sciences Building (LSB) room 213A
(off the lounge) at 2:00-2:30pm Friday September 7th, 2018.

## Background assumed for students
For this class we are not assuming students have background in programming/scripting, nor in bioinformatics. We do assume that students have a working knowledge of basic molecular biology and genetics and have basic familiarity with using computers. i.e. you can figure out how to install basic software on a Mac (OS X) or a PC (Windows).

# What students will need
A laptop with internet access and the ability to install several programs (in particular, `R`, `Python` and a shell emulator in not using a Mac (OS X) or linux.

## Course goals
The primary goal of this course it to provide graduate students an opportunity to develop fundamental computational skills necessary to go on and (in the future) develop the appropriate (and more advanced) skills for bioinformatics, genomics, etc.

## **What this course is not**
Because of limitations of time (one two hour lecture a week for 13-14 weeks), we are purposefully making this a course about fundamental skills. As such, this course will not cover in any detail:

    - Genomic analysis pipelines (RNAseq, variant calling and populations genomics). These are covered in the winter-spring in Bio722.

    - Theory of computer science (nor theory on programming, algorithms, data structures etc).

    - Despite using `R` for much of this course, it is most definitely not a statistics course. Bio708 (taught by Dr. Ben Bolker and Dr. Jonathan Dushoff) is such a course (also using R as the primary programming environmental for statistical modeling.)

    - A bioinformatics class (i.e. we will not teach any conceptual or theoretical background in bioinformatics. All examples will be real examples, but mostly to illustrate the computational skills necessary to run an analysis, not the why).


## Learning Objectives

### Topics (some TBD)
 It is important to note that in order to keep things flexible depending on how things go with the class, these topics are subject to change if necessary. We will discuss in class.

A. Introduction to UNIX and the command line. (Brian)
  1. Introduction to basic shell commands, logging onto remote systems
  2. Standard UNIX utilities that make your day to day computer work (and bioinformatics) easier.
  3. using pipes in UNIX (and the model of streaming data), batch processing of data.
  4. Writing shell tools.
  5. Using your UNIX skills for practical bioinformatic problems (probably setting up a BLAST database, and querying some sequences)
  6. (maybe) Regular expressions are you friend. No really. Using `grep` and its variants (i.e. `agrep`) and `sed` and `awk` for file manipulation and processing

B. Fundamentals of programming using `R`(Ian). [Link to R portal for class](https://github.com/DworkinLab/Bio720/blob/master/Introduction_to_R.md)
  1. Fundamentals of programming in R.
  2. How to avoid repetitive strain injury while programming. Control flow in R (`for` loops, `if else`, etc). Using the `apply` family of functions in `R`. Simple  simulations.
  3. Working with data in `R`. Getting data in. Data munging (subsetting, merging, cleaning). Working with strings in `R`.
  4. Basics of plotting in `R`. Other topics TBD.
  5. Reproducible research using markdown for reports and git for version control.
  6. An Introduction to bioinformatic tools in R. Primarily an introduction to BioConductor, and genomic range data.

C. **This will likely not be taught this year** Fundamentals of program using `python`.

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
You are responsible for ordering your own copies of these books. Both are excellent with only a small amount of overlap, but we are only highly recommending the first book (BDS) for this class. The reason for this is that this year we are only using `UNIX` (and shell scripting) and programming in `R` which are both covered a bit in the BDS book.

[Bioinformatics Data Skills](http://www.amazon.ca/Bioinformatics-Data-Skills-Reproducible-Research/dp/1449367372/ref=sr_1_1?s=books&ie=UTF8&qid=1440614667&sr=1-1&keywords=bioinformatics+data+skills). **HIGHLY RECOMMENDED** This is a new book, but which fills an important gap in that is oriented towards the day to day skills for anyone working in the fields of genomics and bioinformatics. In addition to covering the basic UNIX skills (and why we use UNIX in bioinformatics and genomics), it also covers subjects like overviews of the essential file types (`.fasta`, `.fastq`, `.gff`, etc) that are ubiquitous in the field. There is also a nice, but brief introduction to the essentials of `R`, using bioconductor and in particular range data, and two important chapters on how to organizing (and maximize reproducibility) of computational projects. Currently (August 29th 2018) this is ~52.44$ on amazon.ca . It is available as an [e-book as well from the publisher](http://shop.oreilly.com/product/0636920030157.do). The author is still a PhD student (in population genomics), and wrote this in their first year of graduate school, so definitely worth supporting.

[Practical Computing for Biologists](http://www.amazon.com/s/ref=nb_sb_ss_c_0_24?url=search-alias%3Dstripbooks&field-keywords=practical+computing+for+biologists&sprefix=practical+computing+for+biologists%2Caps%2C144). This book provides a nice, gentle introduction to the basic computational skills all biologists should have. In particular, with introduction to using the UNIX command line, shell scripting, basic `python` programming, regular expressions, working on remote machines and a few other topics. The book is written to be agnostic with respect to discipline (i.e. it is not a bioinformatics book *per se*), but does a great job of being both very accessible and immediately useful. It seems a bit pricy on Amazon.ca, but look around for used copies (it is 4 years old). If you plan to continue in computational research, this is a fantastic resource.


## Important websites

[For Brian's section](http://helix.mcmaster.ca/720.html). This will have pertinent links to Brian's section of the course.

[R tutorial screencasts](https://github.com/DworkinLab/Bio720/blob/master/Introduction_to_R.md). A link to the playlist for screencasts I have put together for the `R` tutorials. I will be adding more as the semester progresses. Mostly we will be using the excellent[Datacamp](https://www.datacamp.com/) online interactive 'courses' for the introductory stuff, and moving on from there.
