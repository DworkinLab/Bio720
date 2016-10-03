#Introduction to Programming with R

Most of the activities for this section of the course will be self-guided, using a combination of video tutorials and practical exercises. In class we will do  It is expected you will have done all of these for class, as we will be doing a larger in class set of activities assuming a basic level of comfort and familiarity. These will then be followed by independent exercises that will be graded.

**Readings**: Unfortunately Practical Computing for Biologists does not have an introduction to `R`. If you did pick up Bioinformatics Data Skills, chapter 8 presents a crash course in programming in `R`. We will not be using any of the functionality in `ggplot2` quite yet, so you can skip those pages (i.e. skip 207-215, 224-227). However, there are plently of amazing [`R` tutorials](https://cran.r-project.org/other-docs.html) (including less amazing ones by me) online. Use the link in the previous sentence to get to some links. I suggest something like [this one](https://cran.r-project.org/doc/contrib/Paradis-rdebuts_en.pdf) to start with. Alternatively McMaster has access to a number of good ebooks. I think that [this one](http://catalogue.mcmaster.ca/catalogue/Record/2702791) by Larry Pace & Joshua Wiley is a good one to start with. You should (while on campus) be able to access the PDF from [here](http://link.springer.com/book/10.1007%2F978-1-4842-0373-6). I am also making a longer list of various `R` resources to help you out, which is currently [here](https://github.com/DworkinLab/Bio720/blob/master/Some_R_resources.md).  I suggest that reading a bit, and/or going through the relevant video tutorial & exercises, and then moving onto the next component may work best.

## Installing R
If you do not have a fairly recent version of `R` installed on your local computer (V.3.3.1 or newer), this is required to be able to complete the class activities. There are several versions of `R` you might consider. R works on all major operating systems (Windows, Apple OS X and Linux).
- You can [download](http://cran.utstat.utoronto.ca/) and install the version of `R` appropriate to your computer. For Mac OS X or Windows you can download them at the page above. For Linux, use `yum`, `apt` or other package management utility you like. For the Mac OS X R GUI, it has a simple script editor that does syntax highlighting, and displays argument flags for functions. I think the Windows R script editor is much more bare bones.
- Alternatively, you can use [R-studio](https://www.rstudio.com/), which is a pretty nice IDE (integrated development environment) for `R`, including advanced syntax highlighting (including RMarkdown, which we will use), and integration with github for version control. While I have a few pet peaves with it, some definitely prefer it, especially when getting started.

## R video tutorials and exercises.
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


## Notes for Ian to add & improve the screencasts (students can ignore this for now).
- seperate "getting data into R" into 3-4 shorter screencasts(with exercises)
- seperate " control flow" into one on if & ifelse, and another on loops.
- seperate the apply screencasts (one for apply alone, additional one for tapply, sapply, lapply, aggregate)
- new screencast on lists as a class (heterogeneous collections) and working with them
- new screencast on data.frame as a list, but also relationship to matrix.
- new screencast on ordering data sets, and the use of the index in R.
- introduce tidyr, dplyr, reshape2, data.table, readr, libraries
- base plotting Introduction
- ggplot2 plotting syntax (and the Hadleyverse syntax as a "macro language" within R)
- real programming things? (declarative syntax, S3 and S4 objects,  )
- things to make R work better with big data sets (pre-allocating space for incoming data, pre-allocating objects and filling them instead of growing objects. Dynamic typing/checking at runtime)
- knitr (with and without Rstudio)
- version control and github (with and without Rstudio)
