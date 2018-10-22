---
title: "Bio720 Introduction to `R`, in class exercise"
author: "Ian Dworkin"
date: "October 22, 2018"
output: 
  html_document: 
    keep_md: yes
    number_sections: yes
    toc: yes
---
## Overview
In class tonight we are going to both practice some of the `R` skills you were introduced to in data camp. You will be working in pairs.

The learning objectives for today are as follows:

1. Learn some best practices for organizing computational projects (i.e. any projects with some scripts, data and outputs). *We already went over this*
2. A brief introduction to markdown and R-markdown to help make your research (and your in class assignments) reproducible and clear. 
3. Learn a little bit about R as a programming language, and where R fits into the ecosystem of programming languages (what it is best to use it for, and why).
4. Learn some intuitive (but not necessarily technical) ideas about *data structures* in general, and review some of the data structures in `R`. 
5. Practice some of the skills that were introduced in the datacamp tutorials. In particular to get a sense that while there are seemingly many ways to accomplish the same task in `R`, they are not always equivalent with respect to speed (and sometimes what they do to the attributes of objects).
6. Even more practice, time permitting!

## A quick introduction to Markdown (in the How to organize computational projects page).

Please [click here](https://github.com/DworkinLab/Bio720/blob/master/IntroductionMarkdownAndVersionControl/Bio720_IntroductionMarkdown.md#what-is-markdown) to link to a brief discussion on these points.

Question 1 - Create your first markdown document! 

In RStudio, go to File > New file > R Markdown. This will give you a few options (which you can explore on your own. Some useful tricks for making a presentation that is independent of powerpoint or Beamer.). Give it a name (something Pithy like your_initials_Bio720_InClassWeek6), and keep the default output as .html. Look in the folder where you saved all of this. What kinds of files do you see?

Click on this, and go through the code. See how code is embedded in the the triple "```" followed by code and ending with another "```"

While (once you have saved the `.Rmd` file) you can run `knitr()` directly from the console, but we will use the RStudio *graphical user interface* (GUI), click "Knit", and then "Knit to html". You will need to save this in a file (I suggest giving it the same name. Use underscore, no spaces!)


Now let's edit this a bit. You need to compute the minimum and maximum values of the `speed` variable in the `cars` object (use just a single function to do this). Then also compute the standard deviation `sd()` of speed.  Before the call to the function add some text describing what you are trying to do with the code. Make it very clear what you are doing by using emphasis (italics) and bold. Then re-knit the document. 


While the code embedding is generic, you can specify that you want to use code from another language. For instance you can add (this may not work on windows...)


```bash
pwd
```

```
## /Users/ian/TeachingAndLectures/Bio720/Bio720/R_exercises
```

You can also do this via the "Insert" button.

Finally, I want you to change some of the default options for the knitting of the document.  Click the options button (looks like a weird starfish). Click on output options. Check the box to include a table of contents. Hit advanced, and click on "keep markdown source file". Save the .Rmd file. re-knit, and take a look in the folder, any additional files? (textedit or notepad should be fine to open these) What are the differences between the .md and the .Rmd file?

## Where does `R` belong in the ecosystem of programming languages.

At a very simple level most (if not all) computer programming languages are ["Turing Complete"](https://en.wikipedia.org/wiki/Turing_completeness), which for our purposes mean that they can in principle be used to do most programming tasks. However in some languages how any particular task is programmed (and how efficient it does the task) varies considerably. Thus some programming languages are better suited for some tasks VS. others. To get a sense of where `R` fits into this ecosystem (i.e. when and where to use `R`) we need to learn a few things. So time to use your google-fu.

Question 2. Is R a low-level or high level programming language? What is the difference?


Question 3. Is R a compiled or interpreted language? What is the difference?

Question 4. R is often described as "slow" or "inefficient". Why? Why would a language with these attributes be so popular? In other words, what were the major goals for the `R` programming language.

Question 5. What computer programming languages is R primarily built in? Why is it not all in `R`

Question 6. Hadley Wickham (the author of many R packages and the R super guru) points out that R is slow in part because of the R language *definition* (such as it is) and in part because of the most commonly used *implementation* `gnu-R`. In very simple terms or analogy descrive what is meant by the *implementation* of a computer language as compared to the *definition*?

 ** Now you are a computer wizard **  Woot Woot!

## Some very basic thoughts on *data structures* in `R`
We are not going to have a computer science-esque discussion of data structures (there are whole courses on this), but instead try to introduce a few basic concepts to understand why computers need to store different types of data in different ways, and why we need to be aware of that.

### What is the point of data structures? (class discussion)
- What kind of data do we want the computer to store for us?
- Why does it matter what kind of data (integers, floating point numbers, strings, Boolean,...)?

### Data structures in R

As was discussed in the datacamp video screencasts, R has a number of different basic data structures (and more can be made or extended by programmers like you!). We need to start with the so-called *atomic* types that can be stored as vectors. 

Question 7. What are the *atomic* types in `R`?

Let's think about a few of these basic atomic types:



```r
x <- 1
```

You can find out information about this with a variety of functions:


```r
str(x)
```

```
##  num 1
```

```r
mode(x)
```

```
## [1] "numeric"
```

```r
typeof(x)
```

```
## [1] "double"
```
Question 8? Why do `mode(x)` and `typeof(x)` give different results?


Let's create a few more objects and take a look at them


```r
y = c(3, 4, 5)
```

Will `x` and `y` differ?  Check and explain?

Now let's create a new object z:


```r
z = 3:5
```
Question 9: How should `y` and `z` compare? how would you check? How can you compare them to see if `R` is treating them the same?

Question 10. Note the behaviour here between `z` and some mathematical operations with `z`. Explain what is happening and why.


```r
typeof(z)
```

```
## [1] "integer"
```

```r
typeof(z+z) # addition
```

```
## [1] "integer"
```

```r
typeof(z-z) # subtraction
```

```
## [1] "integer"
```

```r
typeof(z/z) # regular division
```

```
## [1] "double"
```

```r
typeof(z*z) # regular multiplication
```

```
## [1] "integer"
```

```r
typeof(z^2) # taking the square
```

```
## [1] "double"
```

```r
typeof(sqrt(z))
```

```
## [1] "double"
```

### Some other atomic types in R
Ok, let's think about some of the other basic data types we learned about (strings or "character" in R, boolean/logical)


First let's clean up our workspace and re-assign our variables.


```r
rm( x, y, z) # clean up
x = 1
y = "1"
z = "one"
a = TRUE
b = "TRUE"
```

Question 11. Before checking, think about what types each of these objects should be? Then check. Which (if any) should be identical to one another? Why? How about `y` and `z`? `a` and `b`?


Question 12: How would you get `R` to coerce `x` and `y` to be exactly the same? How about `a` and `b`? 

### Making variables of a specific type

So if you want to make sure you are generating variables of a specific type, you can use the functions `as.x()` where x is the type of variable. You can also specify them at the beginning.

So if I wanted to specify a vector is of type `double` (which is what numeric is)


```r
double_up <- double(length = 10)
identical(numeric(10), double(10))
```

```
## [1] TRUE
```

Question 13. How would I specify a vector of integers of length 10? Why are all the numbers for this (or the line of code above) 0?

### Boolean/logical
Using TRUE/FALSE (logical, Boolean) as an atomic type is quite important. Often we need to check equality of objects, or elements within our objects (considering vectors are a basic storage type in R). These are particularly (as we will see later) useful when subsetting using the index of a vector (or matrix or data.table).

Let's clean up a bit:


```r
rm(x, y, z, a, b)
```

As we will see, both T or TRUE can be used. Likewise for F or FALSE. It is also important to know you don't want these to be treated as strings/characters so don't put quotes around them.


```r
x = T
x
```

```
## [1] TRUE
```

```r
y = TRUE
y
```

```
## [1] TRUE
```

```r
x == y
```

```
## [1] TRUE
```

```r
identical(x, y)
```

```
## [1] TRUE
```

It is also useful to know that TRUE has a numeric value associated with it (1), and FALSE is associated with 0. 


```r
sum(x)
```

```
## [1] 1
```

```r
as.numeric(x)
```

```
## [1] 1
```

```r
a = F
sum(a)
```

```
## [1] 0
```

Question 14: Before running the code to find out, what will the sum of the following vector be `c(rep(T, 10), rep(F, 18), rep(T, 10), rep(F, 6))`

## Building up our data structures. 

Now that we have some better idea (hopefully) of some of the atomic data types, we want to use these to build more complex data structures that may (eventually - like next week) be useful for data analysis, simulations and the like. There are a few important ones that we will use a lot: matrix, list, data.frame, factors, and formula (which we will not cover in Bio720 but is essential for statistical analyses). There are other important ones (like array) but we will cover these other ones first.

Before we get any further and create some new objects, how do we see all of the objects we currently have in our global environment?


```r
ls()
```

```
## [1] "a"         "double_up" "x"         "y"
```

Let's work with a clean slate. How might we remove all of the objects and start fresh? Obviously you could just do a `rm()` command with each object name, but you can also remove all at once.


```r
rm(list=ls())
ls()
```

```
## character(0)
```

Q15. Describe what this command has done. 


Now we are going to create a few new objects and use these to examine some of the properties of our more complex data structures


```r
gene1 <- c(3, 4, 7, 9, 12, 6)
gene2 <- c(11, 17, 12, 25, 23, 7)
gene3 <- c(100, 103, 97, 94, 106, 111)
```
What mode and type should these objects be?


## understanding `factors` in R.

Question 15. Create an object `genotype` of length 6, where the first three observations have the value "wildtype" and the last three are "mutant"


Compare these different objects, genotype (or genotype2 which is identical) and genotype3 (using gl). Are they the same?


So let's think about what a factor is?

factors in R may appear as `character` but for efficiency are stored as integers. The idea is you will have far fewer factor levels (which you can check with `nlevels()`) than number of observations, so this can save memory and speed up computation. However, this means you need to realize that factors are not a special form of `character`, but a special form of `numeric`!

Question 16. If we wanted to make genotype2 into a factor (we will call it genotype2_factor) how would you do so? Is this the same as making it a factor from the very beginning?

How about if we wanted to make genotype3 into a character vector?


```r
genotype3_character <- as.character(genotype3)
```

```
## Error in eval(expr, envir, enclos): object 'genotype3' not found
```

```r
genotype3_character 
```

```
## Error in eval(expr, envir, enclos): object 'genotype3_character' not found
```

```r
class(genotype3_character)
```

```
## Error in eval(expr, envir, enclos): object 'genotype3_character' not found
```

```r
mode(genotype3_character)
```

```
## Error in mode(genotype3_character): object 'genotype3_character' not found
```

```r
identical(genotype3_character, genotype2)
```

```
## Error in identical(genotype3_character, genotype2): object 'genotype3_character' not found
```

```r
genotype3_character == genotype2
```

```
## Error in eval(expr, envir, enclos): object 'genotype3_character' not found
```

Question 17. Let's say we had a second experimental factor which was the day of sampling (3,6) but we want to treat it as a factor `c(3, 6, 3, 6, 3, 6)` how would you code this?

Question 18. What happens if you coerce `day` into a character? Seemingly strange behaviour? However think about it for a minute and try to explain it.


```r
as.character(day)
```

```
## Error in eval(expr, envir, enclos): object 'day' not found
```

Question 19. How about if you coerce day into numeric?

```r
as.numeric(day)
```

```
## Error in eval(expr, envir, enclos): object 'day' not found
```


Question 20. So if you want to turn these into the numbers 3 and 6, how would you do it?

## Back to our data structures of interest. 

Question 21. Provide two different ways of combining `gene1`, `gene2` and `gene3` into a matrix (gene_mat1 and gene_mat2)? Are these the same?

Question 22. How might you fix the issue we observed?

Question 23. Let's take our (character) vectors for day and genotype and use `cbind()` (treatments). Before starting write down whether you think the object `treatments` will have class `matrix`. What will the mode be? Why?

Question 24. Now let's take all of  objects that are vectors of different atomic types (gene1, gene2, gene3, genotype, day) and use `cbind` on them. Call this object `all_the_data`. Before writing the code, write down what you think the class of the object will be. How about the mode/type of the elements of `all_the_data`?

Question 25. Explain why `all_the_data` is the class and has the mode that it does?

## data structures with heterogeneous objects.
 Clearly we did not want to produce a matrix of strings. So we need some sort of data structures where elements (at least at the level of individual vectors that are being organized together) can be of different atomic types (i.e. a collection of heterogeneous objects). There are two main approaches to this, one is the data.frame, which is the spreadsheet like object that is often the easiest to work with when importing data (for analysis and plotting). THe other is a list. As I mentioned in the video tutorials, the data.frame is really a special kind of list. However it is worth comparing and contrasting both. First remove the old `all_the_data` object and make a new one that is a data frame.

## `data.frames`
First let's make a data.frame: 


```r
rm(all_the_data)
```

```
## Warning in rm(all_the_data): object 'all_the_data' not found
```

```r
all_the_data <- data.frame(gene1, gene2, gene3, genotype, day)
```

```
## Error in data.frame(gene1, gene2, gene3, genotype, day): object 'genotype' not found
```

```r
str(all_the_data)
```

```
## Error in str(all_the_data): object 'all_the_data' not found
```

```r
class(all_the_data)
```

```
## Error in eval(expr, envir, enclos): object 'all_the_data' not found
```

```r
mode(all_the_data)
```

```
## Error in mode(all_the_data): object 'all_the_data' not found
```

What class is `all_the_data`? How about `mode`? What is going on?

Notice a couple of interesting thing. First it's class is a data.frame, but it is actually a list underneath. Second, without asking or warning us, it has coerced *genotype* and *day* into factors. It is assuming that since you are treating this like regular data (that you will probably want to analyze or plot) you want these as factors. Often this is true. If you don't want this behaviour there is an argument that you can set `stringsAsFactors == FALSE`.

It is really important to note that data.frames are useful for heterogeneous objects **ONLY IF** all objects (vectors) are the **same length**. It is ok to have missing data, as long as R knows there should be missing data (NA) in certain spots. When you need to store a collection of heterogeneous objects, but the objects are of different lengths, then you need to use lists.

As we showed in the video tutorials and exercises, you can extract and subset in a couple of ways (like lists or as a matrix). So show three different ways to extract the 2nd, 3rd and 4th column from `all_the_data` 

Question 26. Using standard numeric subsetting, extract columns 2,3 and 4.


Question 27. Second, subset using the names of the columns:


Question 28. Third, by using the extraction operator `$` which is how you extract elements from lists:

This approach is useful for single columns, not so useful when you want to extract a bunch though.




Take home message: This can definitely get confusing, but different programmers still use each, so it is important to recognize what is the same and what is different. Importantly the single `[` can be used to extract more than one element from the object, while `$` and `[[` can only select a single element at a time. 


Hopefully it is pretty clear, but there are a couple of ways of adding on new variables. Let's create a gene4 (also of length 6) and add it to the `all_the_data` data.frame


```r
all_the_data$gene4 <- c(10, 11, 7, 11, 2, 3)
```

```
## Error in all_the_data$gene4 <- c(10, 11, 7, 11, 2, 3): object 'all_the_data' not found
```

```r
all_the_data
```

```
## Error in eval(expr, envir, enclos): object 'all_the_data' not found
```

```r
str(all_the_data)
```

```
## Error in str(all_the_data): object 'all_the_data' not found
```

### Lists in `R`

When you need to store a collection of heterogeneous objects, but the objects are of different lengths, then you need to use lists. As you saw above with the `$` and the `[[` operators, you can extract things from lists. However, making lists is simpler than unmaking (well unlisting) lists as we will see.

Make a list called `list_the_data` using the same objects that were used to make `all_the_data`. What will the class of the object be? how about the mode of the objects within the list?


```r
list_the_data = list(gene1, gene2, gene3, genotype, day)
```

```
## Error in eval(expr, envir, enclos): object 'genotype' not found
```

```r
list_the_data
```

```
## Error in eval(expr, envir, enclos): object 'list_the_data' not found
```

```r
str(list_the_data)
```

```
## Error in str(list_the_data): object 'list_the_data' not found
```

```r
names(list_the_data)
```

```
## Error in eval(expr, envir, enclos): object 'list_the_data' not found
```
A couple of things of note:
- It should be pretty clear that there is something different about the way it is storing the information. The names (which is an attribute) of the original objects have been lost. 
- Also, as a list, it does not make assumptions about how you will use the underlying objects, so it has not coerced the character vectors to factors. 

Question 27. How might we get the names of the underlying objects? 

Annoyingly:

```r
list_the_data = list(gene1 = gene1, gene2 = gene2, gene3 = gene3, genotype = genotype, day = day)
```

```
## Error in eval(expr, envir, enclos): object 'genotype' not found
```

```r
list_the_data
```

```
## Error in eval(expr, envir, enclos): object 'list_the_data' not found
```

```r
str(list_the_data)
```

```
## Error in str(list_the_data): object 'list_the_data' not found
```

```r
names(list_the_data)
```

```
## Error in eval(expr, envir, enclos): object 'list_the_data' not found
```

We can also have objects that have different lengths within the list


```r
list_the_data$random_variable = c(T,T,F) 
```

```
## Error in list_the_data$random_variable = c(T, T, F): object 'list_the_data' not found
```

```r
list_the_data
```

```
## Error in eval(expr, envir, enclos): object 'list_the_data' not found
```

```r
str(list_the_data)
```

```
## Error in str(list_the_data): object 'list_the_data' not found
```
 
We can extract variables from lists in a slight variant of the approach we have used so far

 

```r
list_the_data$gene1
```

```
## Error in eval(expr, envir, enclos): object 'list_the_data' not found
```

```r
list_the_data[1]
```

```
## Error in eval(expr, envir, enclos): object 'list_the_data' not found
```

```r
list_the_data["gene1"]
```

```
## Error in eval(expr, envir, enclos): object 'list_the_data' not found
```

```r
list_the_data[[1]]
```

```
## Error in eval(expr, envir, enclos): object 'list_the_data' not found
```

```r
list_the_data[["gene1"]]
```

```
## Error in eval(expr, envir, enclos): object 'list_the_data' not found
```

However, these objects are not all equivalent


```r
class(list_the_data$gene1)
```

```
## Error in eval(expr, envir, enclos): object 'list_the_data' not found
```

```r
class(list_the_data[1])
```

```
## Error in eval(expr, envir, enclos): object 'list_the_data' not found
```

```r
class(list_the_data["gene1"])
```

```
## Error in eval(expr, envir, enclos): object 'list_the_data' not found
```

```r
class(list_the_data[[1]])
```

```
## Error in eval(expr, envir, enclos): object 'list_the_data' not found
```

```r
class(list_the_data[["gene1"]])
```

```
## Error in eval(expr, envir, enclos): object 'list_the_data' not found
```

```r
str(list_the_data$gene1)
```

```
## Error in str(list_the_data$gene1): object 'list_the_data' not found
```

```r
str(list_the_data[1])
```

```
## Error in str(list_the_data[1]): object 'list_the_data' not found
```

```r
str(list_the_data["gene1"])
```

```
## Error in str(list_the_data["gene1"]): object 'list_the_data' not found
```

```r
str(list_the_data[[1]])
```

```
## Error in str(list_the_data[[1]]): object 'list_the_data' not found
```

```r
str(list_the_data[["gene1"]])
```

```
## Error in str(list_the_data[["gene1"]]): object 'list_the_data' not found
```

So using the `[` operator keeps the information (and variable name) as a list, while the `$` or `[[` operators extract just the elements, but do not keep the name. It can not always be coerced in a sensible way as we have done before.

i.e:


```r
as.numeric(list_the_data[1]) 
```

```
## Error in eval(expr, envir, enclos): object 'list_the_data' not found
```

You can `unlist` the vector though

```r
str(as.numeric(unlist(list_the_data[1])))
```

```
## Error in unlist(list_the_data[1]): object 'list_the_data' not found
```

However, this also strips off the name! So you are best to not use the `[` if you can avoid it when using lists. This is not always possible though, so knowing about unlist is useful!


