---
title: "In class assignment week 2, part 2. A worked example using control flow (for loops, if statements, etc)"
author: "Ian Dworkin"
output: 
  html_document:
    keep_md: yes
    number_sections: yes
    toc: yes
---

# Introduction
Let's do a little exercise integrating some of the things we have learned. Here are some Illumina HiSeq reads for one of our recent projects:


```r
read_1 <- "CGCGCAGTAGGGCACATGCCAGGTGTCCGCCACTTGGTGGGCACACAGCCGATGACGAACGGGCTCCTTGACTATAATCTGACCCGTTTGCGTTTGGGTGACCAGGGAGAACTGGTGCTCCTGC"

read_2 <- "AAAAAGCCAACCGAGAAATCCGCCAAGCCTGGCGACAAGAAGCCAGAGCAGAAGAAGACTGCTGCGGCTCCCGCTGCCGGCAAGAAGGAGGCTGCTCCCTCGGCTGCCAAGCCAGCTGCCGCTG"

read_3  <- "CAGCACGGACTGGGGCTTCTTGCCGGCGAGGACCTTCTTCTTGGCATCCTTGCTCTTGGCCTTGGCGGCCGCGGTCGTCTTTACGGCCGCGGGCTTCTTGGCAGCAGCACCGGCGGTCGCTGGC"
```

Question 1. what species are these sequences from?

Question 2. Put all three of these reads into a single object (a vector).  What class will the vector `reads` be? Check to make sure! How many characters are in each read (and why does `length()` not give you what you want.. try...)


```r
reads <- c(read_1, read_2, read_3)
nchar(reads)
```

```
## [1] 124 124 124
```

Question 3. Say we wanted to print each character (not the full string) from read_1, how do we do this using a for loop? You may wish to look at a function like `strsplit()` to accomplish this (there are other ways.)

`strsplit("nameOfString", split ="CharacterVectorToUseForSplitting", fixed = TRUE)`

You might try something like""

```r
for (i in read_1) {
	print(i)
}
```

```
## [1] "CGCGCAGTAGGGCACATGCCAGGTGTCCGCCACTTGGTGGGCACACAGCCGATGACGAACGGGCTCCTTGACTATAATCTGACCCGTTTGCGTTTGGGTGACCAGGGAGAACTGGTGCTCCTGC"
```

While this prints out the full string, it does not do so character by character. THis does work in some other languages though.


We can make this one string into a vector of strings with each element being a single letter


```r
read_1_split <- strsplit(read_1, split = "", fixed = T) # a list
```

Question 4. What kind of object does this return? How might we make it a character vector again?

Since it is a list, which is annoying we want to coerce this into a character vector so we can work with it more easily.


```r
read_1_char <- as.character(unlist(read_1_split))
mode(read_1_char)
```

```
## [1] "character"
```



```r
for (i in read_1_split){
	print(i)
}
```

```
##   [1] "C" "G" "C" "G" "C" "A" "G" "T" "A" "G" "G" "G" "C" "A" "C" "A" "T"
##  [18] "G" "C" "C" "A" "G" "G" "T" "G" "T" "C" "C" "G" "C" "C" "A" "C" "T"
##  [35] "T" "G" "G" "T" "G" "G" "G" "C" "A" "C" "A" "C" "A" "G" "C" "C" "G"
##  [52] "A" "T" "G" "A" "C" "G" "A" "A" "C" "G" "G" "G" "C" "T" "C" "C" "T"
##  [69] "T" "G" "A" "C" "T" "A" "T" "A" "A" "T" "C" "T" "G" "A" "C" "C" "C"
##  [86] "G" "T" "T" "T" "G" "C" "G" "T" "T" "T" "G" "G" "G" "T" "G" "A" "C"
## [103] "C" "A" "G" "G" "G" "A" "G" "A" "A" "C" "T" "G" "G" "T" "G" "C" "T"
## [120] "C" "C" "T" "G" "C"
```

**VS** looping through each character


```r
for (i in read_1_char){
	print(i)
}
```

```
## [1] "C"
## [1] "G"
## [1] "C"
## [1] "G"
## [1] "C"
## [1] "A"
## [1] "G"
## [1] "T"
## [1] "A"
## [1] "G"
## [1] "G"
## [1] "G"
## [1] "C"
## [1] "A"
## [1] "C"
## [1] "A"
## [1] "T"
## [1] "G"
## [1] "C"
## [1] "C"
## [1] "A"
## [1] "G"
## [1] "G"
## [1] "T"
## [1] "G"
## [1] "T"
## [1] "C"
## [1] "C"
## [1] "G"
## [1] "C"
## [1] "C"
## [1] "A"
## [1] "C"
## [1] "T"
## [1] "T"
## [1] "G"
## [1] "G"
## [1] "T"
## [1] "G"
## [1] "G"
## [1] "G"
## [1] "C"
## [1] "A"
## [1] "C"
## [1] "A"
## [1] "C"
## [1] "A"
## [1] "G"
## [1] "C"
## [1] "C"
## [1] "G"
## [1] "A"
## [1] "T"
## [1] "G"
## [1] "A"
## [1] "C"
## [1] "G"
## [1] "A"
## [1] "A"
## [1] "C"
## [1] "G"
## [1] "G"
## [1] "G"
## [1] "C"
## [1] "T"
## [1] "C"
## [1] "C"
## [1] "T"
## [1] "T"
## [1] "G"
## [1] "A"
## [1] "C"
## [1] "T"
## [1] "A"
## [1] "T"
## [1] "A"
## [1] "A"
## [1] "T"
## [1] "C"
## [1] "T"
## [1] "G"
## [1] "A"
## [1] "C"
## [1] "C"
## [1] "C"
## [1] "G"
## [1] "T"
## [1] "T"
## [1] "T"
## [1] "G"
## [1] "C"
## [1] "G"
## [1] "T"
## [1] "T"
## [1] "T"
## [1] "G"
## [1] "G"
## [1] "G"
## [1] "T"
## [1] "G"
## [1] "A"
## [1] "C"
## [1] "C"
## [1] "A"
## [1] "G"
## [1] "G"
## [1] "G"
## [1] "A"
## [1] "G"
## [1] "A"
## [1] "A"
## [1] "C"
## [1] "T"
## [1] "G"
## [1] "G"
## [1] "T"
## [1] "G"
## [1] "C"
## [1] "T"
## [1] "C"
## [1] "C"
## [1] "T"
## [1] "G"
## [1] "C"
```

Question 5. How about if we wanted the number of occurrences of each base? Or their frequencies? You could write a loop, but I suggest looking at the help for the `table()` function... Also keep in mind that for for most objects `length()` tells you how many elements there are in a vector. For lists use `lengths()` (so you can either do this on a character vector or a list, your choice)

For a list you could do this:


```r
site_freq <- table(read_1_split)/lengths(read_1_split)
print(site_freq, digits = 2)
```

```
## read_1_split
##    A    C    G    T 
## 0.19 0.28 0.32 0.21
```
Note the use of `lengths()` not `length()` for the list above!!!!

For a vector of characters you can use this:


```r
site_freq_alt <- table(read_1_char)/length(read_1_char)
print(site_freq_alt, digits = 2)
```

```
## read_1_char
##    A    C    G    T 
## 0.19 0.28 0.32 0.21
```

Since this is a vector of characters we can use `length()`

Question 6. How would you make this into a nice looking function that can work on either  a list or vectors of characters? (Still just for a single read)

You could write a function like:


```r
BaseFrequencies <- function(x) {
    if (mode(x) == "list") {
    	tab <- table(x)/lengths(x)}
    
    else {
    	tab <- table(x)/length(x)
    }	
   return(tab)
}
```


Question 7. Now how can you modify your approach to do it for an arbitrary numbers of reads? You could use a loop or use one of the apply like functions (which one)?

Now try doing it with sapply across the three reads

It may be worth re-writing this function to be a bit more general.


Question 8. Can you revise your function so that it can handle the input of either a string as a single multicharacter vector, **or** a vector of individual characters **or** a list? Try it out with the vector of three sequence reads (`reads`). 


```r
BaseFrequencies <- function(x) {
    
    # if it is a single string still
    if (length(x) == 1 & mode(x) == "character") {
    	x <- strsplit(x, split = "", fixed = T) 
        x <- as.character(unlist(x))
    }     
    
    if (mode(x) == "list") {
    	tab <- table(x)/lengths(x)}
    
    else {
    	tab <- table(x)/length(x)
    }	
   return(tab)
}
```
Let's try it out:


```r
basefreq <- sapply(reads, BaseFrequencies, USE.NAMES = F)
print(basefreq, digits = 2)
```

```
##   [,1]  [,2]  [,3]
## A 0.19 0.274 0.081
## C 0.28 0.331 0.331
## G 0.32 0.298 0.347
## T 0.21 0.097 0.242
```
