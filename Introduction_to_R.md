# Introduction to Programming with R

Most of the activities for this section of the course will be self-guided, using a combination of video tutorials and practical exercises. In class we will do  It is expected you will have done all of these for class, as we will be doing a larger in class set of activities assuming a basic level of comfort and familiarity. These will then be followed by independent exercises that will be graded.

**DataCamp tutorials**: We will be using these a lot in the class! 

You will get email notifications to do the tutorials.



**Readings**: In Bioinformatics Data Skills, chapter 8 presents a crash course in programming in `R`, but this may not be enough detail. We will not be using any of the functionality in `ggplot2` quite yet, so you can skip those pages (i.e. skip 207-215, 224-227). One book which you can access freely is [R for data science](http://r4ds.had.co.nz/index.html) by Garrett Grolemund and Hadley Wickham. It is important to note that it focuses less on programming *per se*, and more on working with data (cleaning, parsing, modifying etc..). It also focuses on using Tidyverse tools in `R` (as opposed to base `R`).

There are plenty of other amazing [`R` tutorials](https://cran.r-project.org/other-docs.html) (including less amazing ones by me) online. Use the link in the previous sentence to get to some links. I suggest something like [this one](https://cran.r-project.org/doc/contrib/Paradis-rdebuts_en.pdf) to start with. Alternatively McMaster has access to a number of good ebooks. I think that [this one](http://catalogue.mcmaster.ca/catalogue/Record/2702791) by Larry Pace & Joshua Wiley is a good one to start with. You should (while on campus) be able to access the PDF from [here](http://link.springer.com/book/10.1007%2F978-1-4842-0373-6). I am also making a longer list of various `R` resources to help you out, which is currently [here](https://github.com/DworkinLab/Bio720/blob/master/Some_R_resources.md).  I suggest that reading a bit, and/or going through the relevant video tutorial & exercises, and then moving onto the next component may work best.

## Installing R

**Even if you have R on your computer, you need to install the current version (and delete your old one). Otherwise not all R libraries will work for you!**

Please install the current version of `R` (4.4.1 "Race for your life" released on 2024/06/14). This is required to be able to complete the class activities. There are several versions of `R` you might consider. R works on all major operating systems (Windows, Apple OS X and Linux). 
- You can [download](https://cran.r-project.org/) and install the version of `R` appropriate to your computer. For Mac OS X or Windows you can download them at the page above. For Linux, use `yum`, `apt` or other package management utility you like. For the Mac OS X R GUI, it has a simple script editor that does syntax highlighting, and displays argument flags for functions. I think the Windows R script editor is much more bare bones.
- You will also need to install the *current* version of [R-studio](https://posit.co/download/rstudio-desktop/), current version (Version: 2024.04.2-764), which is a pretty nice IDE (integrated development environment) for `R`, including advanced syntax highlighting (including RMarkdown, which we will use), and integration with github for version control. While I have a few pet peaves with it, some definitely prefer it, especially when getting started. Certainly makes your research easily reproducible very easily.
- We will be using Git and Github for version control and ultimately assignments in the class. This means you will need to [install git](https://git-scm.com/). I would probably look at [this tutorial](http://happygitwithr.com/install-git.html).
- You will also need to sign up for a free [Github account](https://github.com/). Don't worry about it for now, but you can contact them for an academic account to get private repositories if you want.

## Access to DataCamp tutorials and exercises
I have set up class access to [DataCamp](https://app.datacamp.com/groups/introduction-to-basic-bioinformatics-skills/dashboard) for 2024. You should get an email about this soon, and I will provide details about which courses I expect you to do prior to each class. These will provide foundational skills for us to build on during in class activities, and your problem sets.

## Assignments

**PLEASE NOTE FOR 2024** 
For the second week of class (September 9th, 2024), please complete both the `introduction to R` and the `Intermediate R` courses on DataCamp. 

## Old In class activities (2022)

**For 2024, the R activities are on the main README page**



2.  In class R_exercises for control flow in `R` (as well as the apply family of functions) are  [here](./R_exercises/Bio720_R_InClass_Control_flow_worked_example.Rmd), with [answers](./R_exercises/Bio720_R_InClass_Control_flow_worked_answers.md) and the second activity is [here](./R_exercises/Bio720_R_InClassExercise2.Rmd), with answers [here](./R_exercises/Bio720_R_InClassExercise2Answers.md).

3.  In class exercise on data munging (using base `R`, not with the Wickham-verse libraries) is [here](./R_exercises/Bio720_R_week3_DataMungingInClass.Rmd) and answers [here](./R_exercises/Bio720_R_week3_DataMunging.md).

4.  In class exercises on the basics of plotting (using base `R` not `ggplot2` or `lattice`). [Basic ideas](./R_exercises/Bio720_R_PlottingBasicsInClass.Rmd), and the full code [here](./R_exercises/Bio720_R_PlottingBasics.md) The anatomy of a more advanced example is [here](./R_exercises/AdvancedBasePlotting_CI_Bands.md).
We also introduced some basic concepts about simulations, and went through a simple deterministic example which is [here](./R_exercises/Bio720_SimulatingData_Part1.md). Note that the equations will not display, so download either the [PDF](./R_exercises/Bio720_SimulatingData_Part1.pdf) or [html](./R_exercises/Bio720_SimulatingData_Part1.html) version.  You may also find [this book](https://link.springer.com/book/10.1007%2F978-0-387-89882-7) helpful. McMaster has free access to it.

5. Deterministic and stochastic simulations. markdown is [here](./R_exercises/Bio720_SimulatingData.md), but the [PDF](./R_exercises/Bio720_SimulatingData.pdf) may be easier for the equations. Or download the[.Rmd](./R_exercises/Bio720_SimulatingData.Rmd) and run it yourself.

6. String exercise [here](./R_exercises/Bio720_AnightofStringsInClass.Rmd) and answers [here](./R_exercises/Bio720_AnightofStrings.Rmd)

## other R video tutorials and exercises

Much of this below is from previous versions of the class that belong to the dustbins of (teaching history). However, for posterity (and I don't like getting rid of files), I am keeping them here. There is exactly no need to go through these, unless you have some strange desire to see how I taught some of this more than a decade ago.... 


The data set that is used for some of these activities can be found on the [DRYAD Digital repository](http://datadryad.org/) right [here](http://datadryad.org/bitstream/handle/10255/dryad.8377/dll.csv?sequence=1). You can also set this up (so you do not need a local copy of the data by putting this command in your script or copying and pasting it into the R editor :
```R
dll.data <- read.csv("http://datadryad.org/bitstream/handle/10255/dryad.8377/dll.csv", h=T)
```

*Please note, the scripts may look a bit different, as I have edited them a bit after making the screencasts. All of the important parts are still there!*

The first link is to the screencast itself (hosted on youtube). The subsequent links are to the scripts and exercises. For the first week you only need to watch (and do the exercises) for the first seven parts (up through and including regular sequences and indexing). I am going to make some changes for subsequent weeks.

1. Why use `R` (and why learn to program): Motivating example of working counts of expression data from RNAseq (not yet completed, so skip ahead!).
2. Introduction to `R`: [part 1. `R` as a calculator](https://youtu.be/Kyxx9_NLlUY)
  - the introductory script can be found [here](./Rscripts/R_Introductory_tutorial_part_1.R)
  
3. Introduction to `R`: [part 2. Basic operations and operators in `R`](https://www.youtube.com/watch?v=UrtWeRPpWCw)
 
4. Introduction to `R`: [part 3. Element-by-element operations, booleans & basic functions](https://www.youtube.com/watch?v=8VcysxMmpg0)
  
5. Introduction to `R`: [part 4. objects and classes in `R`](https://www.youtube.com/watch?v=-qDiqnEVaLk)
  
6. Introduction to `R`: [part 5. workspaces and getting help in `R`](https://www.youtube.com/watch?v=0Y9IRfJwzjo)
  
7. Introduction to `R`: [part 6. writing your own functions in `R`](https://www.youtube.com/watch?v=Mth_tvrxik0)
  
8. Introduction to `R`: [part 7. regular sequences and indexing](https://www.youtube.com/watch?v=V5_vb7gLtrk)
  
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
