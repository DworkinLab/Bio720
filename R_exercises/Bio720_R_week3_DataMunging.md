---
title: "Bio720_DataMunging"
author: "Ian Dworkin"
date: "November 5th, 2018"
output: 
  html_document: 
    keep_md: yes
    number_sections: yes
    toc: yes
---


# Data parsing and munging in Base R

## Brief introduction
One of the most common "programming" activities you will do in genomics and bioinfomatics is to get the data into a format that facilitates (and makes possible) the analyses you are doing. This activity can often occupy very large amounts of your time

These activities have many names: *wrangling data, munging data, cleaning data, data manipulation*. Given how important this is, it will come as no surprise that there are entire R libraries devoted to it. There are several lists that are worth looking at like [here](http://datascienceplus.com/best-packages-for-data-manipulation-in-r/) and [here](https://www.analyticsvidhya.com/blog/2015/12/faster-data-manipulation-7-packages/). Also, books like [Data Manipulation with R](https://www.amazon.ca/Data-Manipulation-R-Phil-Spector/dp/0387747303)!

## Why no tidyverse in this tutorial
You have already been introduced to a little bit of the use of the *dplyr* library and verbs within it like `filter()`, `arrange()`, `group_by()` and `summarize()` (one other important one is `select()`). You will have more opportunities to learn aspects of the tidyverse (more on dplyr and ggplot2 in particular). However, today I am not going to show you data munging this way, instead relying on base R techniques.

Why? tidyverse libraries are excellent in many ways, *but* they do have some issues. First and foremost, there are lots of dependencies (so if you load library "a", you also need libraries "b", "c", "d", etc). That can add a lot of libraries potentially. Second, unlike base R where the code is very (probably too) conservative, and those that maintain and alter the code do so in a manner that is very unlikely to alter downstream behaviour (i.e. code in other libraries), tidyverse libraries change more rapidly (and thus, while still very unlikely are more likely to break existing code). 

Currently, *[Bioconductor](http://bioconductor.org/)* objects don't always "play nice" with tidyverse tibbles and approaches. These are essentially two independently evolving "ecosystems of syntax" (really macro languages) in *R*. There is a move to integrate these better. But at this time, if you are planning on using libraries from *Bioconductor* for genomic and bioinformatic analysis it is important to know how to work with data in base R.

It is worth thinking about when to use tidyverse, and when not to.

If you are writing scripts that are just for your own data analysis (i.e. you are not writing a general purpose library) and are not interacting with *Bioconductor* objects too much (or don't mind altering objects from one to the other), then the *tidyverse* is an extremely smart choice. Fast to learn, relatively consistent (in comparison to Base *R*), and often faster than base *R*.

## *data.table* library for fast import and execution of large data sets.
Finally it is worth pointing out that the *data.table* library is not only the most efficient for importing large data sets (by a wide margin usually), but data munging (sorting, subsetting, etc..) can also be done via this library with its own syntax. For large datasets, this is usually the fastest. Datacamp has a course on this [here](https://www.datacamp.com/courses/data-table-data-manipulation-r-tutorial). Like `readr` and *tibbles* in tidyverse, the `data.table` is a slightly different object from a standard `data.frame`. If you are working with large datasets in `R`, knowing this library is essential! 

## for your own review

We have already spent considerable time in previous tutorials making data frames using `cbind`, `rbind`, `data.frame`, etc.. As we have limited time in class I leave review of that to you from the screencasts and activites at datacamp.

## install packages.
You may wish to install certain `R` libraries that get used a lot. While R-studio has a number of easy ways to install libraries I still prefer having it in the script. You only need to install the libraries once, although it is worth updating them regularly, and I re-install them everytime I update `R`.


```r
# remove hash to install.
#install.packages("data.table")
```

## Data we are using.

We will use a number of different data sets to illustrate each of the points, but for the main one, we will use an old *Drosophila melanogaster* data set measuring several traits (lengths) and the number of sex comb teeth (a structure used to clasp females during copulation) for different wild type strains (line) reared at different developmental temperatures, with and without a mutation that effects proximal-distal axis development in limbs (genotype). Note that we are downloading the data from a data repository (Dryad). If you really care, it is from [this](https://www.ncbi.nlm.nih.gov/pubmed/15733306) paper.



```r
dll_data = read.csv("http://datadryad.org/bitstream/handle/10255/dryad.8377/dll.csv", header=TRUE)
```


Before we go on, how should we look at the data to make sure it imported correctly, and the structure (and other information) about the object we have just created?


```r
summary(dll_data)
```

```
##    replicate         line      genotype        temp          femur     
##  Min.   :1.00   line-7 : 132   Dll: 871   Min.   :25.0   Min.   :0.21  
##  1st Qu.:1.00   line-18: 121   wt :1102   1st Qu.:25.0   1st Qu.:0.53  
##  Median :1.00   line-4 : 112              Median :25.0   Median :0.55  
##  Mean   :1.18   line-8 : 110              Mean   :27.4   Mean   :0.55  
##  3rd Qu.:1.00   line-2 : 104              3rd Qu.:30.0   3rd Qu.:0.57  
##  Max.   :2.00   line-11: 100              Max.   :30.0   Max.   :0.70  
##                 (Other):1294                             NA's   :24    
##      tibia          tarsus          SCT      
##  Min.   :0.34   Min.   :0.11   Min.   : 6.0  
##  1st Qu.:0.46   1st Qu.:0.18   1st Qu.:10.0  
##  Median :0.48   Median :0.19   Median :11.0  
##  Mean   :0.48   Mean   :0.19   Mean   :11.2  
##  3rd Qu.:0.50   3rd Qu.:0.20   3rd Qu.:12.0  
##  Max.   :0.61   Max.   :0.26   Max.   :32.0  
##  NA's   :19     NA's   :17     NA's   :25
```

```r
str(dll_data)
```

```
## 'data.frame':	1973 obs. of  8 variables:
##  $ replicate: int  1 1 1 1 1 1 1 1 1 1 ...
##  $ line     : Factor w/ 27 levels "line-1","line-11",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ genotype : Factor w/ 2 levels "Dll","wt": 1 1 1 1 1 1 1 1 1 1 ...
##  $ temp     : int  25 25 25 25 25 25 25 25 25 25 ...
##  $ femur    : num  0.59 0.55 0.588 0.588 0.596 ...
##  $ tibia    : num  0.499 0.501 0.488 0.515 0.502 ...
##  $ tarsus   : num  0.219 0.214 0.211 0.211 0.207 ...
##  $ SCT      : int  9 13 11 NA 12 14 11 12 10 12 ...
```

```r
dim(dll_data)
```

```
## [1] 1973    8
```

```r
head(dll_data)
```

```
##   replicate   line genotype temp femur tibia tarsus SCT
## 1         1 line-1      Dll   25 0.590 0.499  0.219   9
## 2         1 line-1      Dll   25 0.550 0.501  0.214  13
## 3         1 line-1      Dll   25 0.588 0.488  0.211  11
## 4         1 line-1      Dll   25 0.588 0.515  0.211  NA
## 5         1 line-1      Dll   25 0.596 0.502  0.207  12
## 6         1 line-1      Dll   25 0.577 0.499  0.207  14
```
- Since we are using base *R*, this is a regular data.frame.

## Cleaning data
### removing missing data

Sometimes your data set has missing data, i.e. for some reason you could not measure one of your variables on a particular object. How you decide to deal with missing data can be a big topic, but for the moment we are going to assume you want to delete rows that contain missing data. 

First let's check if there is any missing data


```r
head(is.na(dll_data))
```

```
##      replicate  line genotype  temp femur tibia tarsus   SCT
## [1,]     FALSE FALSE    FALSE FALSE FALSE FALSE  FALSE FALSE
## [2,]     FALSE FALSE    FALSE FALSE FALSE FALSE  FALSE FALSE
## [3,]     FALSE FALSE    FALSE FALSE FALSE FALSE  FALSE FALSE
## [4,]     FALSE FALSE    FALSE FALSE FALSE FALSE  FALSE  TRUE
## [5,]     FALSE FALSE    FALSE FALSE FALSE FALSE  FALSE FALSE
## [6,]     FALSE FALSE    FALSE FALSE FALSE FALSE  FALSE FALSE
```

```r
sum(is.na(dll_data))
```

```
## [1] 85
```
So there are rows that contain missing data somewhere. 


More generally we can ask:

```r
anyNA(dll_data)
```

```
## [1] TRUE
```

### Why missing data can cause headaches

```r
mean(dll_data$femur)
```

```
## [1] NA
```
- Since there is some missing data, the `mean()` function defaults to NA, to warn you. This is a useful behaviour.


```r
mean(dll_data$femur, na.rm = TRUE)
```

```
## [1] 0.546
```


### What to do with this rows with missing data
We need to decide what to do with the missing data. If you have lots of data, and do not mind removing all cases (rows) that contain **any** missing data then `na.omit` is a useful solution.


```r
dll_data_complete <- na.omit(dll_data)

dim(dll_data)
```

```
## [1] 1973    8
```

```r
dim(dll_data_complete)
```

```
## [1] 1918    8
```

This has caused us to remove 55 observations. It may be useful to instead look at where the missing data is contained and maybe make a subsetted data set based on that. There are numerous options about how `R` should handle operations if an `NA` is found. You can look at the help under `?NA`. 


For the moment, we are going to accept the approach above, and use it. However, I don't generally recommend this. Instead, you can first subset the variables you are going to use, and then decide how to deal with the missing data. Or you can allow them into your analysis, but make sure to set the right flags in functions that enable them to deal with missing data.


```r
dll_data <- na.omit(dll_data)
mean(dll_data$femur)
```

```
## [1] 0.546
```
- Note that the mean now is not the same as before. This is because additional variables (not just femur) had missing data for different observations, but we removed any observation with any missing data!

### finding and removing duplicate rows in a data frame.
The `duplicated()` allows you to see which, if any rows are perfect duplicates of one another.

We can first ask if any rows are duplicates.

This provides a Boolean

```r
head(duplicated(dll_data))
```

```
## [1] FALSE FALSE FALSE FALSE FALSE FALSE
```

```r
tail(duplicated(dll_data))
```

```
## [1] FALSE FALSE FALSE FALSE FALSE FALSE
```

So we can do something like:


```r
sum(duplicated(dll_data))
```

```
## [1] 0
```

Or if you want it as a Boolean

```r
any(duplicated(dll_data))
```

```
## [1] FALSE
```

`R` has a convenience function for the above:


```r
anyDuplicated(dll_data)
```

```
## [1] 0
```
So for this example there was one duplicate rows.


```r
dll_data[anyDuplicated(dll_data),]
```

```
## [1] replicate line      genotype  temp      femur     tibia     tarsus   
## [8] SCT      
## <0 rows> (or 0-length row.names)
```
Which in this case, is because there is so much missing data for this observation, and only one variable (SCT) had actually data for a discrete value. So it is likely real.

However let's say an accident happened in generating the data set and a few rows were duplicated. I am going to use the `sample()` function to extract 5 random rows.


```r
new_rows <- dll_data[sample(nrow(dll_data), size = 5, replace = T ),]

dll_data2 <- rbind(dll_data, new_rows)
str(dll_data2)
```

```
## 'data.frame':	1923 obs. of  8 variables:
##  $ replicate: int  1 1 1 1 1 1 1 1 1 1 ...
##  $ line     : Factor w/ 27 levels "line-1","line-11",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ genotype : Factor w/ 2 levels "Dll","wt": 1 1 1 1 1 1 1 1 1 1 ...
##  $ temp     : int  25 25 25 25 25 25 25 25 25 25 ...
##  $ femur    : num  0.59 0.55 0.588 0.596 0.577 ...
##  $ tibia    : num  0.499 0.501 0.488 0.502 0.499 ...
##  $ tarsus   : num  0.219 0.214 0.211 0.207 0.207 ...
##  $ SCT      : int  9 13 11 12 14 11 12 10 12 13 ...
##  - attr(*, "na.action")= 'omit' Named int  4 61 73 92 93 142 207 268 315 319 ...
##   ..- attr(*, "names")= chr  "4" "61" "73" "92" ...
```

Now you don't know which rows, they are.  Please show how you would find whether there are any duplicated rows, and which ones they are using `duplicated()`

```r
any(duplicated(dll_data2))
```

```
## [1] TRUE
```

```r
dll_data2[duplicated(dll_data2),]
```

```
##       replicate      line genotype temp femur tibia tarsus SCT
## 14351         1 line-OreR      Dll   30 0.559 0.474  0.166   8
## 11781         1   line-18      Dll   30 0.496 0.440  0.165   8
## 4831          1   line-11       wt   25 0.518 0.484  0.193  12
## 8271          2    line-3       wt   25 0.552 0.471  0.190  10
## 12281         1    line-2      Dll   30 0.490 0.423  0.166   8
```

So what should we do? We can use the `unique()` to remove these duplicated rows if we know for certain that they do not belong in the data set (i.e. they are duplicated in error).


```r
dll_data_unique <- unique(dll_data2)

dim(dll_data_unique)
```

```
## [1] 1918    8
```

```r
dim(dll_data2)
```

```
## [1] 1923    8
```

```r
dim(dll_data)
```

```
## [1] 1918    8
```

Let's just do a bit of clean-up.


```r
rm(dll_data_complete, dll_data_unique, dll_data2)
```

## Subsetting data sets (review)

In some sense this is a review of what we learned during the first week, but there are a few details that are worth considering.

Let's look at the data again:


```r
str(dll_data)
```

```
## 'data.frame':	1918 obs. of  8 variables:
##  $ replicate: int  1 1 1 1 1 1 1 1 1 1 ...
##  $ line     : Factor w/ 27 levels "line-1","line-11",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ genotype : Factor w/ 2 levels "Dll","wt": 1 1 1 1 1 1 1 1 1 1 ...
##  $ temp     : int  25 25 25 25 25 25 25 25 25 25 ...
##  $ femur    : num  0.59 0.55 0.588 0.596 0.577 ...
##  $ tibia    : num  0.499 0.501 0.488 0.502 0.499 ...
##  $ tarsus   : num  0.219 0.214 0.211 0.207 0.207 ...
##  $ SCT      : int  9 13 11 12 14 11 12 10 12 13 ...
##  - attr(*, "na.action")= 'omit' Named int  4 61 73 92 93 142 207 268 315 319 ...
##   ..- attr(*, "names")= chr  "4" "61" "73" "92" ...
```

As we can see there are two genotypes. Let's say we wanted to make a data frame (`dll_data_wt`) that was a subset with only the wild type (`wt`) data? How would you do this with the index?

### using the index

```r
dll_data_wt <- dll_data[dll_data$genotype == "wt",]
```

Has this created the appropriate subset of data?


```r
with(dll_data, table(genotype))
```

```
## genotype
##  Dll   wt 
##  841 1077
```

```r
nrow(dll_data_wt)
```

```
## [1] 1077
```

So it seems like it is working.

### subsetting on factors and removing unused levels.

So that seems to match. We should only have one level now in `genotype`. Let's check.


```r
levels(dll_data_wt$genotype)
```

```
## [1] "Dll" "wt"
```

So what is going on? 

We know we only have the correct number of observations associated, but it still has the level `Dll` in the factor genotype, why? Because of the way it has been stored, each of these levels is still associated as part of `genotype`. Thankfully it is very easy to fix. Use the `droplevels()` to drop unused levels.


```r
dll_data_wt <- droplevels(dll_data_wt)
levels(dll_data_wt$genotype)
```

```
## [1] "wt"
```

### using the subset function

Of course the easiest alternative is to use the `subset` function we introduced a few weeks back. Please use this to make a data frame (`dll_data_Dll`) for the *Dll* factor level of genotype and check what has happened to the *wt* factor level. 



```r
dll_data_Dll <- subset(dll_data, genotype == "Dll" )
dim(dll_data_Dll)
```

```
## [1] 841   8
```

```r
with(dll_data, table(genotype))
```

```
## genotype
##  Dll   wt 
##  841 1077
```

```r
levels(dll_data_Dll$genotype)
```

```
## [1] "Dll" "wt"
```

```r
dll_data_Dll <- droplevels(dll_data_Dll)
```

`subset()` like using the index allows you to select specific columns as well. Create a new version of dll_data_Dll where the only columns that remain are line, genotype, temp and SCT.


```r
dll_data_Dll <- subset(dll_data, 
                       genotype == "Dll",
                       select = c(line, genotype, temp, SCT))

dim(dll_data_Dll)
```

```
## [1] 841   4
```

# the `%in%` matching operator 

`%in%` can be very useful for finding things in your data set to work on.

The `%in%` is a matching operator (also see `match()`) that returns a boolean (TRUE or FALSE). This can be really useful to find things. Go ahead and find all instances where a couple of the lines `c("line-Sam", "line-1")` are in the data set `dll_data$line`. I recommend first creating a dummy logical variable that can be used to filter via the index`[]`.


```r
matched_set <- dll_data$line %in% c("line-Sam", "line-1")
sum(matched_set)
```

```
## [1] 122
```

```r
dll_data_new_subset <- dll_data[matched_set,]
dim(dll_data_new_subset)
```

```
## [1] 122   8
```


### Clean up before moving on.
let's remove the data frames we are not using.


```r
rm(dll_data_Dll, dll_data_wt, dll_data_new_subset, matched_set)
```

## Cleaning variable names
Let's take a closer look at the names of the different fly strains used.


```r
levels(dll_data$line)
```

```
##  [1] "line-1"    "line-11"   "line-12"   "line-13"   "line-15"  
##  [6] "line-16"   "line-17"   "line-18"   "line-19"   "line-2"   
## [11] "line-20"   "line-21"   "line-22"   "line-23"   "line-24"  
## [16] "line-26"   "line-27"   "line-3"    "line-4"    "line-6"   
## [21] "line-7"    "line-8"    "line-9"    "line-CanS" "line-OreR"
## [26] "line-Sam"  "line-w"
```

Some are numbers, some are letters. We decide, that we need to clean this up a bit. First off, we have no need to have "line-" as a prefix for each of these labels, and may make things confusing down the road. So we want to get rid of them. However we **don't** want to edit the spreadsheet, as it is most importantly bad scientific practice, and would be a lot of work. So how do we do it efficiently?

We actually have a couple of options.

First thing to keep in mind though is that `dll_data$line` is currently stored as an object of class `factor`. Why is this important?

### `substr()`
Let's try to do some simple string manipulation. There are a few functions that could be useful here. Let's start with `substr()` (substring). This will extract or replace substrings within a character vector. So what do we need to do first?


```r
line_str <- as.character(dll_data$line)
str(line_str)
```

```
##  chr [1:1918] "line-1" "line-1" "line-1" "line-1" "line-1" "line-1" ...
```

```r
head(line_str)
```

```
## [1] "line-1" "line-1" "line-1" "line-1" "line-1" "line-1"
```

Now we can use `substr()`. It is by position.


```r
line_names <- substr(line_str, 
                     start = 6, stop= 1000000L )

head(line_names)
```

```
## [1] "1" "1" "1" "1" "1" "1"
```

```r
tail(line_names)
```

```
## [1] "w" "w" "w" "w" "w" "w"
```

In this case, since the strings we want to keep start at different positions in the original string, we can't use an arbitrary string length (i.e. try `stop  = 7`). 


### `strsplit()`
An alternative way to do this is to use the fact that we (purposefully) used a delimiter character to seperate parts of the name. In this case we used the "-" to seperate the word "line" from the actual line name we car about. As such we can use a second function that can be very valuable, `strsplit()`. Use `strsplit` to create a new object splitting on the "-"


```r
line_names2 <- strsplit(line_str, split = "-")
head(line_names2)
```

```
## [[1]]
## [1] "line" "1"   
## 
## [[2]]
## [1] "line" "1"   
## 
## [[3]]
## [1] "line" "1"   
## 
## [[4]]
## [1] "line" "1"   
## 
## [[5]]
## [1] "line" "1"   
## 
## [[6]]
## [1] "line" "1"
```
Here it has created a list, where each list contains the 2 elements ("line" and the name of the line which we want). We are going to convert this list to a matrix so we can extract what we need.


```r
line_names2_mat <- matrix(unlist(line_names2), 
                          ncol = 2,
                          byrow = TRUE)

head(line_names2_mat)
```

```
##      [,1]   [,2]
## [1,] "line" "1" 
## [2,] "line" "1" 
## [3,] "line" "1" 
## [4,] "line" "1" 
## [5,] "line" "1" 
## [6,] "line" "1"
```

And the second column contains what we want (line names) so we could make a new variable.


```r
dll_data$line_names <- factor(line_names2_mat[,2])
str(dll_data$line_names)
```

```
##  Factor w/ 27 levels "1","11","12",..: 1 1 1 1 1 1 1 1 1 1 ...
```

As an aside, in my lab I enforce an approach for imaging where all images, videos, or sequence files have a consistent (and identical within project) naming convention. For instance our wing images have a naming convention like initials_project_line_sex_genotype_rep. So for images it may be something like:

`ID_GBE_SAM_F_WT_01.tif`

Which would be ID = Ian Dworkin (initials), GBE = project name, SAM is line name, etc... This enables me to import the file names as variables and use `strsplit()` like we did above and automatically populate the variables for the experiment. No need to try to seperate things in a spreadsheet program!

### `gsub`

Many times however, there is no consistent delimiter like "-", nor can we assume the variable name we want is in a particular part of the string. So sometimes we need to utilize some simple regular expressions. `R` has a number of useful regular expression tools in the `base` library (`gsub`, `grep`, `grepl`, etc ). I suggest using `?grep` to take a deeper look. Here we will encounter just a few.

`sub()` and `gsub()` are functions that perform replacements for the first match or globally (g) for all matches respectively. If we go back to the names of the original lines, we will see they all have `line-` in common, so perhaps we can search and replace based on that pattern. This is quite easy as well. use `gsub` to extract just the line name without the first part of the string, "line-"


```r
str(line_str)
```

```
##  chr [1:1918] "line-1" "line-1" "line-1" "line-1" "line-1" "line-1" ...
```

```r
line_names3 <- gsub(pattern = "line-",
                    replacement = "",
                    x = line_str)

head(line_names3)
```

```
## [1] "1" "1" "1" "1" "1" "1"
```

```r
tail(line_names3)
```

```
## [1] "w" "w" "w" "w" "w" "w"
```

What has happened here? Well we looked in the `line_str` object, and found all instances of our pattern (which in this case was the simple "line-"), and then we *replaced* this with "", in other words an emptry string.

`R` can handle POSIX style regular expressions (probably what you learned with Brian) as well PERL style regular expressions by setting the flag `PERL = TRUE`. This enables very flexible pattern searching.

### Some other things to help "clean" up names, etc..

You may have noticed that for the line names which are short words, some are lower case and some are upper case. It could be that this inconsistentcy make cause an issue down the road. So let's say we wanted to take `line_names3` and make all of the words lower case. How would we do that? `R` has some nice functions to help, namely `tolower()` (and `toupper()`) and the more general `chartr()` for character translation.


```r
line_names3 <- tolower(line_names3)

head(line_names3)
```

```
## [1] "1" "1" "1" "1" "1" "1"
```

```r
tail(line_names3, n = 25)
```

```
##  [1] "sam" "sam" "sam" "sam" "sam" "sam" "w"   "w"   "w"   "w"   "w"  
## [12] "w"   "w"   "w"   "w"   "w"   "w"   "w"   "w"   "w"   "w"   "w"  
## [23] "w"   "w"   "w"
```

Say we wanted to call the "sam" line by its proper name "Samarkand". How could we do this?


```r
line_names3 <- sub("sam", "SAMARKAND", line_names3)
```

#### Clean up before moving on.
rm(line_names2, line_names3, line_names, line_names2_mat, line_str, new_rows)

### renaming variables, factors etc..

For this experiment, all flies were reared during their development at two temperatures, 25C and 30C. Currently this is kept as an integer. 


```r
str(dll_data$temp)
```

```
##  int [1:1918] 25 25 25 25 25 25 25 25 25 25 ...
```

However we might want to treat this as a factor. Also your PI has decided they want you to encode it at "HighTemp" and "LowTemp". So you have decided to create a new variable to do so. Not surprisingly, you have a few different ways to do this. The most obvious (and what you have done before) is to generate a factor. You can make the variable part of the data frame object, but I am going to keep it seperate for now.

Use the `factor()` function to generate a new variable `temp_as_factor` with factor labels "LowTemp", "HighTemp".


```r
temp_as_factor <- with(dll_data,
                      factor(temp, 
                             labels = c("LowTemp", "HighTemp")))

str(temp_as_factor)
```

```
##  Factor w/ 2 levels "LowTemp","HighTemp": 1 1 1 1 1 1 1 1 1 1 ...
```

```r
head(temp_as_factor)
```

```
## [1] LowTemp LowTemp LowTemp LowTemp LowTemp LowTemp
## Levels: LowTemp HighTemp
```

```r
tail(temp_as_factor)
```

```
## [1] HighTemp HighTemp HighTemp HighTemp HighTemp HighTemp
## Levels: LowTemp HighTemp
```

### Using ifelse()

An alternative approach would be to use the `ifelse()` function. Create a new variable (that achieves the same thing) using `ifelse()`


```r
temp_as_factor2 <- with(dll_data, ifelse(temp == 25, "LowTemp", 
                                         ifelse(temp == 30, "HighTemp", NA)))
temp_as_factor2 <- factor(temp_as_factor2)
str(temp_as_factor2)
```

```
##  Factor w/ 2 levels "HighTemp","LowTemp": 2 2 2 2 2 2 2 2 2 2 ...
```

One issue with this is nested `if` statements, like `for` loops can get messy. Plus there is an internal limit on how much nesting can occur (I think a hierarchy of 5 or 6), so you want to take some care in using this approach.

## Sorting data frames

Sometimes you just want to sort your data frame. While this is rarely necessary for statistical analyses or plotting of the data, it is important to look at the raw data as a double (triple!!!) check that there are no typos, or issues you were not aware of. I would say in genomic analysis this is most important when you are looking at the summary of statistical analyses and you want to sort your output by genes with greatest effects or something else. While there is a `sort()` function in `R` this is not generally the way you want to approach it.


```r
head(sort(dll_data$SCT))
```

```
## [1] 6 7 7 7 7 7
```

```r
tail(sort(dll_data$SCT))
```

```
## [1] 18 18 18 18 20 22
```
This just sorts the vector, which does not help us much. Instead we are going to use the `order()` function which provides information on the index for each observation. So 


```r
head(order(dll_data$SCT))
```

```
## [1] 1375 1148 1268 1269 1325 1555
```

```r
tail(order(dll_data$SCT))
```

```
## [1] 1068 1283 1284 1336   60 1430
```
Provides the row numbers. We can use this to sort via the *index*. Use the `order()` to sort the dll_data frame into a new object (dll_data_sorted) based on dll_data$SCT.


```r
dll_data_sorted <- dll_data[order(dll_data$SCT),]

head(dll_data_sorted)
```

```
##      replicate    line genotype temp femur tibia tarsus SCT line_names
## 1409         1  line-8      Dll   30 0.512 0.468  0.169   6          8
## 1179         1 line-18      Dll   30 0.527 0.435  0.164   7         18
## 1300         1 line-23      Dll   30 0.508 0.454  0.170   7         23
## 1301         1 line-23      Dll   30 0.484 0.422  0.168   7         23
## 1358         1  line-4      Dll   30 0.541 0.470  0.158   7          4
## 1601         1 line-16       wt   30 0.559 0.508  0.179   7         16
```

How would you sort on tibia?

### sorting from highest to lowest

So this by default sorts from lowest to highest. Use the help function for `order` and sort the data from highest to lowest.


```r
dll_data_sorted <- dll_data[order(-dll_data$SCT),]
head(dll_data_sorted)
```

```
##      replicate    line genotype temp femur tibia tarsus SCT line_names
## 1469         1  line-w      Dll   30 0.506 0.451  0.174  22          w
## 62           1 line-13      Dll   25 0.583 0.526  0.217  20         13
## 1099         2 line-13      Dll   30 0.589 0.535  0.169  18         13
## 1315         1 line-24      Dll   30 0.540 0.470  0.191  18         24
## 1316         1 line-24      Dll   30 0.534 0.480  0.187  18         24
## 1369         1  line-6      Dll   30 0.504 0.509  0.175  18          6
```

### Sorting on multiple different columns. 

This is a natural extension of what we have already done. Sort first on *SCT*, then on *temp*


```r
dll_data_sorted <- dll_data[order(dll_data$SCT, dll_data$temp),]
head(dll_data_sorted)
```

```
##      replicate    line genotype temp femur tibia tarsus SCT line_names
## 1409         1  line-8      Dll   30 0.512 0.468  0.169   6          8
## 1179         1 line-18      Dll   30 0.527 0.435  0.164   7         18
## 1300         1 line-23      Dll   30 0.508 0.454  0.170   7         23
## 1301         1 line-23      Dll   30 0.484 0.422  0.168   7         23
## 1358         1  line-4      Dll   30 0.541 0.470  0.158   7          4
## 1601         1 line-16       wt   30 0.559 0.508  0.179   7         16
```

```r
head(dll_data)
```

```
##   replicate   line genotype temp femur tibia tarsus SCT line_names
## 1         1 line-1      Dll   25 0.590 0.499  0.219   9          1
## 2         1 line-1      Dll   25 0.550 0.501  0.214  13          1
## 3         1 line-1      Dll   25 0.588 0.488  0.211  11          1
## 5         1 line-1      Dll   25 0.596 0.502  0.207  12          1
## 6         1 line-1      Dll   25 0.577 0.499  0.207  14          1
## 7         1 line-1      Dll   25 0.618 0.494  0.204  11          1
```

## Merging data frames

We often end up in a situation where we collect more data on the same set of samples (or in the same locales) that we need to include in our analyses. Again, we do not want to go back and edit the spreadsheets, but instead merge the data frames during the analysis. For instance, let's say we had this new data that related the elevation that each of the 27 lines were originally collected at, and we thought this may explain some important component of the variation in the traits (in this case size). So our collaborator collected this data:


```r
line_names <- dll_data$line
levels(line_names)
```

```
##  [1] "line-1"    "line-11"   "line-12"   "line-13"   "line-15"  
##  [6] "line-16"   "line-17"   "line-18"   "line-19"   "line-2"   
## [11] "line-20"   "line-21"   "line-22"   "line-23"   "line-24"  
## [16] "line-26"   "line-27"   "line-3"    "line-4"    "line-6"   
## [21] "line-7"    "line-8"    "line-9"    "line-CanS" "line-OreR"
## [26] "line-Sam"  "line-w"
```

```r
elevations <- c(100, 300, 270, 250, 500, 900, 500, 1100, 500,
                3000,500, 570, 150, 800, 600, 500, 1900, 100,
                300, 270, 250, 500, 900, 500, 1100, 500, 600)

MeanDayTimeTemp <- c(rnorm(27, mean = 20, sd = 5))
elevation_data <- data.frame(levels(line_names), 
                             elevations,
                             MeanDayTimeTemp)
```
So how do we combine these?

Using a for loop or ifelse is possible, but would be messy. There are two easier options. First using merge.

### `merge()`

The `merge()` function is designed for exactly this. It works on matching variable names (exactly is best) from variables with the same name, or where the matching pairs are specified.

In our case


```r
names(dll_data)
```

```
## [1] "replicate"  "line"       "genotype"   "temp"       "femur"     
## [6] "tibia"      "tarsus"     "SCT"        "line_names"
```

```r
names(elevation_data)
```

```
## [1] "levels.line_names." "elevations"         "MeanDayTimeTemp"
```

The easiest thing to do would be to rename the variable in `elevation_data` so that it matches that in `dll_data`


```r
names(elevation_data)[1] <- "line"
str(elevation_data)
```

```
## 'data.frame':	27 obs. of  3 variables:
##  $ line           : Factor w/ 27 levels "line-1","line-11",..: 1 2 3 4 5 6 7 8 9 10 ...
##  $ elevations     : num  100 300 270 250 500 900 500 1100 500 3000 ...
##  $ MeanDayTimeTemp: num  22.4 15.7 24.1 17.5 29.6 ...
```

Now we can merge the data


```r
merged_data <- merge(x = elevation_data, 
                     y = dll_data)
```

How can we check if it merged correctly?


```r
head(merged_data)
```

```
##     line elevations MeanDayTimeTemp replicate genotype temp femur tibia
## 1 line-1        100            22.4         1      Dll   25 0.590 0.499
## 2 line-1        100            22.4         1      Dll   25 0.550 0.501
## 3 line-1        100            22.4         1      Dll   25 0.588 0.488
## 4 line-1        100            22.4         1      Dll   25 0.596 0.502
## 5 line-1        100            22.4         1      Dll   25 0.577 0.499
## 6 line-1        100            22.4         1      Dll   25 0.618 0.494
##   tarsus SCT line_names
## 1  0.219   9          1
## 2  0.214  13          1
## 3  0.211  11          1
## 4  0.207  12          1
## 5  0.207  14          1
## 6  0.204  11          1
```

```r
tail(merged_data)
```

```
##        line elevations MeanDayTimeTemp replicate genotype temp femur tibia
## 1913 line-w        600            21.5         1       wt   30 0.544 0.495
## 1914 line-w        600            21.5         1       wt   30 0.507 0.480
## 1915 line-w        600            21.5         1       wt   30 0.528 0.485
## 1916 line-w        600            21.5         1       wt   30 0.545 0.486
## 1917 line-w        600            21.5         1       wt   30 0.542 0.504
## 1918 line-w        600            21.5         1       wt   30 0.511 0.469
##      tarsus SCT line_names
## 1913  0.181  12          w
## 1914  0.181  10          w
## 1915  0.178   9          w
## 1916  0.174  12          w
## 1917  0.167  12          w
## 1918  0.165  11          w
```

The `merge()` function has a great deal more flexibility, and it is worth looking at the examples in the help function! A word of warning with `merge()` and `reshape()` (below). The functions can result in some odd behaviour when there are non-unique rows which create conflicts. They can also stumble with missing data sometimes. So be on the look out for such things when you are using these functions.

### the power of the index

We can also use the index to do something similar.

Say we measured gene expression of a particular transcript for both genotypes *Dll* and *wt*. We wanted to include that. It would be very easy to use an ifelse, but we can also use the index


```r
geneXpression <- c(Dll=1121, wt = 2051)
merged_data$geneExp <- geneXpression[dll_data$genotype]

head(dll_data)
```

```
##   replicate   line genotype temp femur tibia tarsus SCT line_names
## 1         1 line-1      Dll   25 0.590 0.499  0.219   9          1
## 2         1 line-1      Dll   25 0.550 0.501  0.214  13          1
## 3         1 line-1      Dll   25 0.588 0.488  0.211  11          1
## 5         1 line-1      Dll   25 0.596 0.502  0.207  12          1
## 6         1 line-1      Dll   25 0.577 0.499  0.207  14          1
## 7         1 line-1      Dll   25 0.618 0.494  0.204  11          1
```

Make a new variable based on genotype for `geneY = c(200, 500)`.

#### Clean up for next part.

```r
rm(elevation_data, elevations, MeanDayTimeTemp, temp_as_factor, dll_data_sorted, geneXpression, temp_as_factor2, line_names, stuff, merged_data)
```

```
## Warning in rm(elevation_data, elevations, MeanDayTimeTemp,
## temp_as_factor, : object 'stuff' not found
```

## reshaping data

There are often times when your data is in a format that is not appropriate for the analyses you need to do. In particular sometimes you need each variable (say each gene whose expression you are monitoring) in its own column (so called *wide* format). Sometimes it is best to have a single column for gene expression, with a second column indexing which gene the measure of expression is for. i.e. all expression data is in a single column, and thus each row has a single measure, but the number of rows will be equal to the number of samples multiplied by the number of genes (lots of rows!!). This is called the *long* format.

This is such a common operation that `R` has a function `reshape()` to do this. Unfortunately (but perhaps well deserved), the arguments and argument names for the `reshape()` function are not so intuitive. Thus while I give an example of it below, this is one of those cases that using an additional package (like `reshape2` or `tidyr`) may be a worthwhile investment in time. Indeed `gather` and `spread` are really straightforward in *tidyr*

### `reshape()` function in `R`
Let's take a quick look at our data again:

```r
head(dll_data)
```

```
##   replicate   line genotype temp femur tibia tarsus SCT line_names
## 1         1 line-1      Dll   25 0.590 0.499  0.219   9          1
## 2         1 line-1      Dll   25 0.550 0.501  0.214  13          1
## 3         1 line-1      Dll   25 0.588 0.488  0.211  11          1
## 5         1 line-1      Dll   25 0.596 0.502  0.207  12          1
## 6         1 line-1      Dll   25 0.577 0.499  0.207  14          1
## 7         1 line-1      Dll   25 0.618 0.494  0.204  11          1
```

It is clear that we have three variables (femur, tibia and tarsus) which are all measures of size (3 segments of the leg). While we currently have these in the *wide* format, for some analysis (such as in mixed models) we may wish to put them in a *long* format. 

First let me show you the code (and that it worked), and then we will go through the various arguments to make sense of what we have done.

```r
long_data <- reshape(dll_data,
                     varying = list(names(dll_data)[5:7]),
                     direction = "long",
                     ids = row.names(dll_data),
                     times = c("femur", "tibia", "tarsus"),
                     timevar = "leg_segments",
                     v.names = "length")
```

Let's check that this worked. How might we do this? What do we expect in terms of number of rows in `long_data`?


```r
# It should have 3 times the number of rows as our original data
(3*nrow(dll_data)) == nrow(long_data)
```

```
## [1] TRUE
```

```r
head(long_data, n = 10)
```

```
##          replicate   line genotype temp SCT line_names leg_segments length
## 1.femur          1 line-1      Dll   25   9          1        femur  0.590
## 2.femur          1 line-1      Dll   25  13          1        femur  0.550
## 3.femur          1 line-1      Dll   25  11          1        femur  0.588
## 5.femur          1 line-1      Dll   25  12          1        femur  0.596
## 6.femur          1 line-1      Dll   25  14          1        femur  0.577
## 7.femur          1 line-1      Dll   25  11          1        femur  0.618
## 8.femur          1 line-1      Dll   25  12          1        femur  0.582
## 9.femur          1 line-1      Dll   25  10          1        femur  0.572
## 10.femur         1 line-1      Dll   25  12          1        femur  0.568
## 11.femur         1 line-1      Dll   25  13          1        femur  0.555
##          id
## 1.femur   1
## 2.femur   2
## 3.femur   3
## 5.femur   5
## 6.femur   6
## 7.femur   7
## 8.femur   8
## 9.femur   9
## 10.femur 10
## 11.femur 11
```

```r
tail(long_data, n = 10)
```

```
##             replicate   line genotype temp SCT line_names leg_segments
## 1963.tarsus         1 line-w       wt   30  11          w       tarsus
## 1964.tarsus         1 line-w       wt   30  10          w       tarsus
## 1965.tarsus         1 line-w       wt   30  10          w       tarsus
## 1966.tarsus         1 line-w       wt   30  10          w       tarsus
## 1967.tarsus         1 line-w       wt   30  12          w       tarsus
## 1968.tarsus         1 line-w       wt   30  10          w       tarsus
## 1969.tarsus         1 line-w       wt   30   9          w       tarsus
## 1971.tarsus         1 line-w       wt   30  12          w       tarsus
## 1972.tarsus         1 line-w       wt   30  12          w       tarsus
## 1973.tarsus         1 line-w       wt   30  11          w       tarsus
##             length   id
## 1963.tarsus  0.185 1963
## 1964.tarsus  0.185 1964
## 1965.tarsus  0.184 1965
## 1966.tarsus  0.183 1966
## 1967.tarsus  0.181 1967
## 1968.tarsus  0.181 1968
## 1969.tarsus  0.178 1969
## 1971.tarsus  0.174 1971
## 1972.tarsus  0.167 1972
## 1973.tarsus  0.165 1973
```



### some other approaches to using `reshape()`

Sometimes we actually have unique subject identifiers already in the data frame, and we don't need to just use row numbers. These unique identifiers can be a combination of several variables if it generates a unique combination for each row (for data currently in wide format). See the section below on combining names. In turns out that even if we combined line, genotype, temp we would not have unique names, so that strategy does not work here. So we can make a fake identifer

```r
dll_data$subject <- as.character(1:nrow(dll_data))
```


```r
long_data2 <- reshape(dll_data,
                     varying = list(names(dll_data)[5:7]),
                     direction = "long",
                     idvar = "subject",
                     #ids = row.names(dll_data),
                     times = c("femur", "tibia", "tarsus"),
                     timevar = "leg_segments",
                     v.names = "length"
                     )
```

Which we can check again


```r
3*nrow(dll_data) == nrow(long_data2)
```

```
## [1] TRUE
```

```r
head(long_data2)
```

```
##         replicate   line genotype temp SCT line_names subject leg_segments
## 1.femur         1 line-1      Dll   25   9          1       1        femur
## 2.femur         1 line-1      Dll   25  13          1       2        femur
## 3.femur         1 line-1      Dll   25  11          1       3        femur
## 4.femur         1 line-1      Dll   25  12          1       4        femur
## 5.femur         1 line-1      Dll   25  14          1       5        femur
## 6.femur         1 line-1      Dll   25  11          1       6        femur
##         length
## 1.femur  0.590
## 2.femur  0.550
## 3.femur  0.588
## 4.femur  0.596
## 5.femur  0.577
## 6.femur  0.618
```

In this case, this is really similar to the approach above as we used the same fundamental IDs.

## Combining variables together.
Sometimes you need to combine several variables together. For instance, we may need to look at combinations of both *genotype* and *temp*


```r
str(dll_data)
```

```
## 'data.frame':	1918 obs. of  10 variables:
##  $ replicate : int  1 1 1 1 1 1 1 1 1 1 ...
##  $ line      : Factor w/ 27 levels "line-1","line-11",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ genotype  : Factor w/ 2 levels "Dll","wt": 1 1 1 1 1 1 1 1 1 1 ...
##  $ temp      : int  25 25 25 25 25 25 25 25 25 25 ...
##  $ femur     : num  0.59 0.55 0.588 0.596 0.577 ...
##  $ tibia     : num  0.499 0.501 0.488 0.502 0.499 ...
##  $ tarsus    : num  0.219 0.214 0.211 0.207 0.207 ...
##  $ SCT       : int  9 13 11 12 14 11 12 10 12 13 ...
##  $ line_names: Factor w/ 27 levels "1","11","12",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ subject   : chr  "1" "2" "3" "4" ...
##  - attr(*, "na.action")= 'omit' Named int  4 61 73 92 93 142 207 268 315 319 ...
##   ..- attr(*, "names")= chr  "4" "61" "73" "92" ...
```

What if we wanted to combine these two variables together into a meta-variable? There are a few good options. One is to use paste of course


```r
meta_variable <- with(dll_data,
                      paste(genotype, temp, 
                            sep = ":"))

head(meta_variable, n = 20)
```

```
##  [1] "Dll:25" "Dll:25" "Dll:25" "Dll:25" "Dll:25" "Dll:25" "Dll:25"
##  [8] "Dll:25" "Dll:25" "Dll:25" "Dll:25" "Dll:25" "Dll:25" "Dll:25"
## [15] "Dll:25" "Dll:25" "Dll:25" "Dll:25" "Dll:25" "Dll:25"
```

```r
str(meta_variable)
```

```
##  chr [1:1918] "Dll:25" "Dll:25" "Dll:25" "Dll:25" "Dll:25" "Dll:25" ...
```
Note: you get to choose the character to use as a seperator. ":", "." and "_" are very common!

The only problem is this is now a character vector not a factor, although this is easily fixed.


```r
meta_variable <- as.factor(meta_variable)
str(meta_variable)
```

```
##  Factor w/ 4 levels "Dll:25","Dll:30",..: 1 1 1 1 1 1 1 1 1 1 ...
```
Great!

The other pretty easy way, especially for factors is to use the `interaction` function.


```r
meta_variable2 <- with(dll_data,
     interaction(as.factor(temp), genotype,
                 drop = TRUE,
                 sep = ":"))

head(meta_variable2, n = 20)
```

```
##  [1] 25:Dll 25:Dll 25:Dll 25:Dll 25:Dll 25:Dll 25:Dll 25:Dll 25:Dll 25:Dll
## [11] 25:Dll 25:Dll 25:Dll 25:Dll 25:Dll 25:Dll 25:Dll 25:Dll 25:Dll 25:Dll
## Levels: 25:Dll 30:Dll 25:wt 30:wt
```

```r
str(meta_variable2)
```

```
##  Factor w/ 4 levels "25:Dll","30:Dll",..: 1 1 1 1 1 1 1 1 1 1 ...
```

Please note the `drop = TRUE` flag. You need this if there are combinations of factors that do not exist in the data, and it is almost always the safer choice unless you need all possible combinations (including those that don't occur).

## Please ignore
In class we are going to go through just a few activities that are important. combining levels in a factor (grepl factor, ), transform, numeric to factor (cut),   creating a dataset from multiple external files (RNAseq example)....
