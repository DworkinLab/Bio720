---
title: "BIO720_2021_FunctionsAssignment"
author: "ID"
date: "16/11/2021"
output:
  html_document:
    toc: yes
    number_sections: yes
    keep_md: yes
  pdf_document:
    toc: yes
editor_options:
  chunk_output_type: console
---






# Functions assignment



**Due Date**: Tuesday November 23rd, 2021 (5PM)

This assignment is based on the class material from November 15th, 2021 (ID)

## DataCamp

1. DataCamp course "Intermediate R" Chapters 1- 4. Send me a picture showing you have completed the appropriate chapters, but I should be able to see the completion status on DataCamp as well.

## R functions to write

**Note:** Feel free to either use base R functionality or tidyverse syntax for any of the questions below. 

2. Write a named function that can take a single numeric vector as an input (argument), and outputs (as a *numeric vector*) the sample size, arithmetic mean, standard deviation, standard error of the mean and the coefficient of variation. Make sure to include in the function ways of dealing with unexpected input values like **a)** missing data in the input vector, **b)** an empty input vector, **c)** non-numeric values in the input vector and **d)** provides a warning if the sample size is too small in the input vector (a sample size of 5 or less). Make sure it computes the *correct sign (+)* for the coefficient of variation. Finally make sure the output of this function is clear (i.e. we know which elements of the vector correspond to the mean, standard deviation, etc).

3. Write a named function that can take two named input numeric vectors as arguments. In the body of the function compute the element by element ratios from the two vectors, and then computes the arithmetic mean of those ratios. This function also should compute the ratio of the means of the two vectors. Then, as a list output an object with the following elements: **a** the mean of the element-by-element ratios, **b** the ratio of the means and **c** the difference between these two values. Given the goals of this function (computing these three elements), use the ideas we discussed in class and what you used in the previous function to build in appropriate "safety checks" and warnings against unexpected inputs, or inputs that would result in erroneous output.

4. Now make a revised version of the function you wrote for question 3, where the user of the function can choose to use either the mean or the median (but the functions default behaviour is the mean like in the function you wrote for question 3).

5.  Read in the sequence read below (from an Illumina HiSeq run). Write a function that can take a sequence read like this as input and outputs a table of the frequencies of A, C, G & T (two digits for each) along with the number of base pairs in the read. Keep in mind a few things. You may wish to split the read from a single string into a vector where each element is a single base pair (or into a list if you prefer). There are numerous options but one approach would be to use the `strsplit` function. This will output a list (not a character vector). If you want to use the length function note that it won't have the expected behaviour on a list (but if you look at the help for length you will see an option for this). Some other functions that *may* help you on your way (use R help to understand how they can be used) include: `table`, `nchar`, `length`, `lengths`, `print`, `unlist`, `mode`.




```r
read_1 <- "CGCGCAGTAGGGCACATGCCAGGTGTCCGCCACTTGGTGGGCACACAGCCGATGACGAACGGGCTCCTTGACTATAATCTGACCCGTTTGCGTTTGGGTGACCAGGGAGAACTGGTGCTCCTGC"
```


6. Can you write a revised version of your function from question 5  so that it can handle the input of either a string as a single multi-character vector (like `read_1`), **or** a character vector of individual nucelotide bases **or** a list with each element being a character for an individual base (like you would get from using `strsplit()` on `read_1`)? The `mode` function may be useful for this question (depending on how you choose to code it).

7. While you don't need to write the actual changes to the function from question 6 (but feel free to do so), provide at least three things you might build checks or warnings for that given the function, and possibly "incorrect" inputs. i.e some things a user of the function may try to use as input that would result in unexpected or incorrect outputs.

