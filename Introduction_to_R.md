# Introduction to Programming with R

Most of the activities for this section of the course will be self-guided, using a combination of video tutorials and practical exercises. In class we will do  It is expected you will have done all of these for class, as we will be doing a larger in class set of activities assuming a basic level of comfort and familiarity. These will then be followed by independent exercises that will be graded.

**DataCamp tutorials**: We will be using these a lot in the class!

**Readings**: In Bioinformatics Data Skills, chapter 8 presents a crash course in programming in `R`, but this may not be enough detail. We will not be using any of the functionality in `ggplot2` quite yet, so you can skip those pages (i.e. skip 207-215, 224-227). One book which you can access freely is [R for data science](http://r4ds.had.co.nz/index.html) by Garrett Grolemund and Hadley Wickham. It is important to note that it focuses less on programming *per se*, and more on working with data (cleaning, parsing, modifying etc..). It also focuses on using Tidyverse tools in `R` (as opposed to base `R`).

There are plenty of other amazing [`R` tutorials](https://cran.r-project.org/other-docs.html) (including less amazing ones by me) online. Use the link in the previous sentence to get to some links. I suggest something like [this one](https://cran.r-project.org/doc/contrib/Paradis-rdebuts_en.pdf) to start with. Alternatively McMaster has access to a number of good ebooks. I think that [this one](http://catalogue.mcmaster.ca/catalogue/Record/2702791) by Larry Pace & Joshua Wiley is a good one to start with. You should (while on campus) be able to access the PDF from [here](http://link.springer.com/book/10.1007%2F978-1-4842-0373-6). I am also making a longer list of various `R` resources to help you out, which is currently [here](https://github.com/DworkinLab/Bio720/blob/master/Some_R_resources.md).  I suggest that reading a bit, and/or going through the relevant video tutorial & exercises, and then moving onto the next component may work best.

## Installing R
Please install the current version of `R` (V.3.5.1 or newer). This is required to be able to complete the class activities. There are several versions of `R` you might consider. R works on all major operating systems (Windows, Apple OS X and Linux).
- You can [download](http://cran.utstat.utoronto.ca/) and install the version of `R` appropriate to your computer. For Mac OS X or Windows you can download them at the page above. For Linux, use `yum`, `apt` or other package management utility you like. For the Mac OS X R GUI, it has a simple script editor that does syntax highlighting, and displays argument flags for functions. I think the Windows R script editor is much more bare bones.
- You will also need to install the *current* version of [R-studio](https://www.rstudio.com/), version 1.1.456, which is a pretty nice IDE (integrated development environment) for `R`, including advanced syntax highlighting (including RMarkdown, which we will use), and integration with github for version control. While I have a few pet peaves with it, some definitely prefer it, especially when getting started. Certainly makes your research easily reproducible very easily.
- We will be using Git and Github for version control and ultimately assignments in the class. This means you will need to [install git](https://git-scm.com/). I would probably look at [this tutorial](http://happygitwithr.com/install-git.html).
- You will also need to sign up for a free [Github account](https://github.com/). Don't worry about it for now, but you can contact them for an academic account to get private repositories if you want.

## Access to DataCamp tutorials and exercises
I have set up class access to [DataCamp](https://www.datacamp.com/groups/bio720) for 2018. You should get an email about this soon, and I will provide details about which courses I expect you to do prior to each class. These will provide foundational skills for us to build on during in class activities, and your problem sets.

## Assignments
[R Assignment 1](./assignments/Bio720_R_Assignment1.md) - Due October 21st and 22nd (see notes)

[R Assignment 2](./assignments/Bio720_R_Assignment2.md) - Due October 28th & 29th 2018 (see notes).

[R Assignment 3](./assignments/Bio720_R_Assignment3.md) - Due November 5th before class. Answers are [here](./assignments/Bio720_R_Assignment3_answers.md)

## In class activities
1. In class R exercises on the introduction to `R` are [here](https://github.com/DworkinLab/Bio720/blob/master/R_exercises/R_ClassExercise_1_InClass.md). The [answers](./R_exercises/R_ClassExercise_1_answers.md) Please go through DataCamp courses (intro to R and intermediate R) screencasts 1 - 7 above first (or else it may not make much sense).

2. (October 29th 2018) In class R_exercises for control flow in `R` (as well as the apply family of functions) are  [here](./R_exercises/Bio720_R_InClass_Control_flow_worked_example.Rmd), with [answers](./R_exercises/Bio720_R_InClass_Control_flow_worked_answers.md) and the second activity is [here](./R_exercises/Bio720_R_InClassExercise2.Rmd), with answers [here](./R_exercises/Bio720_R_InClassExercise2Answers.md).

3. (November 5th 2018) In class exercise on data munging (using base `R`, not with the Wickham-verse libraries) is [here](./R_exercises/Bio720_R_week3_DataMungingInClass.Rmd) and answers [here](./R_exercises/Bio720_R_week3_DataMunging.md).

4. (November 19th 2018) In class exercises on the basics of plotting (using base `R` not `ggplot2` or `lattice`). [Basic ideas](./R_exercises/Bio720_R_PlottingBasicsInClass.Rmd), and the anatomy of a more advanced example is [here](./R_exercises/AdvancedBasePlotting_CI_Bands.md)

## other R video tutorials and exercises.
The data set that is used for some of these activities can be found on the [DRYAD Digital repository](http://datadryad.org/) right [here](http://datadryad.org/bitstream/handle/10255/dryad.8377/dll.csv?sequence=1). You can also set this up (so you do not need a local copy of the data by putting this command in your script or copying and pasting it into the R editor :
```R
dll.data <- read.csv("http://datadryad.org/bitstream/handle/10255/dryad.8377/dll.csv", h=T)
```

*Please note, the scripts may look a bit different, as I have edited them a bit after making the screencasts. All of the important parts are still there!*

The first link is to the screencast itself (hosted on youtube). The subsequent links are to the scripts and exercises. For the first week you only need to watch (and do the exercises) for the first seven parts (up through and including regular sequences and indexing). I am going to make some changes for subsequent weeks.

1. Why use `R` (and why learn to program): Motivating example of working counts of expression data from RNAseq (not yet completed, so skip ahead!).
2. Introduction to `R`: [part 1. `R` as a calculator](https://youtu.be/Kyxx9_NLlUY)
  - the introductory script can be found [here](./Rscripts/R_Introductory_tutorial_part_1.R)
  - the first exercise is [here](./R_exercises/R_exercise_1.md)
3. Introduction to `R`: [part 2. Basic operations and operators in `R`](https://www.youtube.com/watch?v=UrtWeRPpWCw)
  - the second exercise is [here](./R_exercises/R_exercise_2.md)
4. Introduction to `R`: [part 3. Element-by-element operations, booleans & basic functions](https://www.youtube.com/watch?v=8VcysxMmpg0)
  - the third exercise is [here](./R_exercises/R_exercise_3.md)
5. Introduction to `R`: [part 4. objects and classes in `R`](https://www.youtube.com/watch?v=-qDiqnEVaLk)
  - the fourth exercise is [here](./R_exercises/R_exercise_4.md)
6. Introduction to `R`: [part 5. workspaces and getting help in `R`](https://www.youtube.com/watch?v=0Y9IRfJwzjo)
  - the fifth exercise is [here](./R_exercises/R_exercise_5.md)
7. Introduction to `R`: [part 6. writing your own functions in `R`](https://www.youtube.com/watch?v=Mth_tvrxik0)
  - the sixth exercise is [here](./R_exercises/R_exercise_6.md)
8. Introduction to `R`: [part 7. regular sequences and indexing](https://www.youtube.com/watch?v=V5_vb7gLtrk)
  - the seventh exercise is [here](./R_exercises/R_exercise_7.md)
9. Introduction to `R`: [part 8. getting data into `R`](https://www.youtube.com/watch?v=SlupCvzH2nM)
  - the script for this activity can be found [here](./Rscripts/R_Introductory_tutorial_part_2.R)
  - the data can be downloaded [here](http://datadryad.org/bitstream/handle/10255/dryad.8377/dll.csv)
10. Introduction to `R`: [part 9. control flow in `R`](https://www.youtube.com/watch?v=FgtqJ8-DN7k)
  - The script can be viewed [here](./Rscripts/IntroductionControlFlowR.R)
11. Introduction to `R`: [part 10. using the apply family of functions in `R`](https://www.youtube.com/watch?v=uL_LdYS-scQ)
  - The script can be viewed [here](./Rscripts/applyLikeFunctionsR.R)


## Notes for Ian to add & improve the screencasts (students can ignore this).
- separate "getting data into R" into 3-4 shorter screencasts(with exercises)
- separate " control flow" into one on if & ifelse, and another on loops.
- separate the apply screencasts (one for apply alone, additional one for tapply, sapply, lapply, aggregate)
- new screencast on lists as a class (heterogeneous collections) and working with them
- new screencast on data.frame as a list, but also relationship to matrix.
- new screencast on ordering data sets, and the use of the index in R.
- introduce tidyr, dplyr, reshape2, data.table, readr, purrr
- base plotting Introduction
- ggplot2 plotting syntax (and the Hadleyverse syntax as a "macro language" within R)
- real programming things? (declarative syntax, S3 and S4 objects,  )
- things to make R work better with big data sets (pre-allocating space for incoming data, pre-allocating objects and filling them instead of growing objects. Dynamic typing/checking at runtime)
- knitr (with and without Rstudio)
- version control and github (with and without Rstudio)
