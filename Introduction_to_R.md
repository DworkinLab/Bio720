#Introduction to Programming with R

 We do not have enough "in class" time to really develop all of the `R` programming skills necessary to enable you to independently perform analysis of differential expression using RNAseq data. Therefore, most of the activities will be self-guided, using a combination of video tutorials and practical exercises. It is expected you will have done all of these for class, as we will be doing a larger in class set of activities assuming a basic level of comfort and familiarity.

Readings: From Bioinformatics Data Skills, chapter 8 presents a crash course in programming in `R`. We will not be using any of the functionality in `ggplot2` quite yet, so you can skip those pages (i.e. skip 207-215, 224-227). I suggest that reading a bit, and/or going through the relevant video tutorial & exercises, and then moving onto the next component may work best.

## Installing R
If you do not have a fairly recent version of `R` installed on your local computer (V.3.2.1 or newer), this is required to be able to complete the class activities. There are several versions of `R` you might consider.
- You can [download](http://cran.utstat.utoronto.ca/) and install the version of `R` appropriate to your computer. For Mac OS X or Windows you can download them at the page above. For Linux, use `yum`, `apt` or other package management utility you like. For the Mac OS X R GUI, it has a simple script editor that does syntax highlighting, and displays argument flags for functions. I think the Windows R script editor is much more bare bones.
- Alternatively, you can use [R-studio](https://www.rstudio.com/), which is a pretty nice IDE (integrayed development environment) for `R`, including advanced syntax highlighting (including RMarkdown, which we will use), and integration with github for version control. While I have a few pet peaves with it, most folks love it.

## R video tutorials and exercises.
The data set that is used for some of these activities can be found on the [DRYAD Digital repository](http://datadryad.org/) right [here](http://datadryad.org/bitstream/handle/10255/dryad.8377/dll.csv?sequence=1). You can also set this up (so you do not need a local copy of the data by putting this command in your script or copying and pasting it into the R editor :
```R
dll.data <- read.csv("http://datadryad.org/bitstream/handle/10255/dryad.8377/dll.csv", h=T)
```

*Please note, the scripts may look a bit different, as I have edited them a bit after making the screencasts. All of the important parts are still there!*

The first link is to the screencast itself (hosted on youtube). The subsequent links are to the scripts and exercises.

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
9. Introduction to `R`: [part 8. getting data into `R`](https://www.youtube.com/watch?v=SlupCvzH2nM)
  - the script for this activity can be found [here](./Rscripts/R_Introductory_tutorial_part_2.R)
  - the data can be downloaded [here](http://datadryad.org/bitstream/handle/10255/dryad.8377/dll.csv)
10. Introduction to `R`: [part 9. control flow in `R`](https://www.youtube.com/watch?v=FgtqJ8-DN7k)
  - The script can be viewed [here](./Rscripts/IntroductionControlFlowR.R)
11. Introduction to `R`: [part 10. using the apply family of functions in `R`](https://www.youtube.com/watch?v=uL_LdYS-scQ) 
  - The script can be viewed [here](./Rscripts/applyLikeFunctionsR.R)
