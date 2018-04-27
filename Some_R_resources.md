
Here I am just providing some links to various `R` resources that may be helpful as you start to learn to program and perform analyses in `R`.


## R *cheat sheets* (i.e quick summaries of key functions)
Not surprisingly, R studio has some to [check out](https://www.rstudio.com/resources/cheatsheets/)


## R Tutorials
There are many of these, but here are a few...

[DataCamp Tutorials](https://www.datacamp.com/home) has nice very interactive introductory tutorials for `R`. In particular the introduction to R, intermediate R, and functions in R "classes" will get you a long way. It has similarly nice  introductory `Python` tutorials. I have set up a class page for Bio720 (2016) so you can use the tutorials.

[R studio tutorials](https://www.rstudio.com/online-learning/)

[swirl interactive R tutorials](http://swirlstats.com/)

[R tutorials on R tutor](http://www.r-tutor.com/)

### Some of my favourite (general) R books
(less about understanding and performing analyses, more about programming and using R). There are many many more, but these are just a few I have used and read.

Check back on [this page](https://github.com/DworkinLab/Bio720/blob/master/Introduction_to_R.md) for a few suggestions of introductory books.

[R programming for bioinformatics](http://www.crcnetbase.com/isbn/9781420063684). Robert Gentleman. A very nice introduction to programming in R. Not much to do about bioinformatics *per se*, but a good way of learning `R` with bioinformatics examples, and some introduction to some of the classes particular to bioconductor.

[The R Book](http://www.mcmu.eblib.com.libaccess.lib.mcmaster.ca/patron/FullRecord.aspx?p=297479). A really nice introduction to using R, data structures in R, plotting and basic statistics. If you are going to be using R, you want to have a copy of this to use!

[The art of R programming](https://www.nostarch.com/artofr.htm). A really nice introduction to programming in R. Does not assume you have any background in programming.

[R in a nutshell](http://shop.oreilly.com/product/0636920022008.do). Once you have some basic understanding of programming in general, this is a nice place to look things up quickly and does a good job of explaining what R is doing with respect to data structures.

[R for Data Science](http://r4ds.had.co.nz/). An excellent introduction to importing, munging and plotting data the **tidyverse** way.

## R syntax style guides
[R syntax style guide according to google](https://google.github.io/styleguide/Rguide.xml)

[R syntax style guide according to Hadley Wickham (and thus R-studio)](http://adv-r.had.co.nz/Style.html)

[R syntax style guide according to me](https://msu.edu/~idworkin/ZOL851_style_guide.html)


## how to find libraries that you might need?
[CHECK OUT CRAN task views](https://cran.r-project.org/web/views/) to see links to libraries with specialized functions.

### For genomics and bioinformatics most R libraries
are available through [bioconductor](http://bioconductor.org/). There are hundreds of amazing libraries. SUper useful tutorials online and example workflows will get you going quickly. One thing to know is that some of the way that data is stored for some bioconductor objects may seem a bit different than other R objects. In particular the two main data structures in bioconductor are the `ExpressionSet` and `SummarizedExperiments`.  It takes a bit to get used to it, but they can be very helpful down the road!


### The *tidyverse*
*The Hadley Wickham-verse* (also known as the [tidyverse](http://tidyverse.org/)). [Hadley Wickham](http://hadley.nz/) has developed some very important and widely used libraries in R, and is one of the chief scientific officers for R-studio. In particular one of his goals was to make doing routine things (data munging, aggregating, cleaning and plotting) in R simpler, more intuitive and more consistent than is found with standard R syntax, and as such many of his libraries are widely used (https://github.com/hadley). Much of the usage of this is clearly demonstrated in the [R for data science](http://r4ds.had.co.nz/) book.

 It is worth noting that of his (HW) effort to do so, he has created a sort of macro language in R (more like a parallel syntax really) that can be quite different. However, while it can be quite different than base R syntax, the learning curve is far less steep, and far more consistent in usage (part of the point). These packages include: dplyr, reshape2, tidyr, readr, devtools, stringr, ggplot2, lubridate, [broom](https://github.com/tidyverse/broom). All of these are available from his [github](https://github.com/hadley/) page. He has also written several books that cover the gamut of library specific things (ggplot2) to advanced R programming and writing your own packages.

 - While the links above should get you started, everything is also at the github page [here](https://github.com/tidyverse). There also seems to be some other libraries that seem to link closely to the tidyverse philosophy like [`tidytext`](https://www.rdocumentation.org/packages/tidytext/versions/0.1.2) for text mining.
data.table has some amazingly useful functions to import and work with large data.

- One thing to note is that the tidyverse and bioconductor represent two fairly distinct ecosystems in `R`, see this [blog post](https://4dpiecharts.com/2017/04/20/mabbles-the-missing-piece-of-a-unified-theory-of-tidyness/) for some points of comparison (and how to create some tidyverse ideas for bioconductor).
Also check out [this list](https://support.rstudio.com/hc/en-us/articles/201057987-Quick-list-of-useful-R-packages) of other suggested packages (some are the same). There is also the [`biobroom`](https://github.com/dgrtwo/biobroom) library for "tidying" bioconductor objects.

### other important libraries.
libraries like [`data.table`](https://github.com/Rdatatable/data.table/wiki) provides a very easy to use approach to reading in large data sets (common in genomics) **very quickly**. `readr` uses a different approach but with similar goals.

## Making reports
In addition to great tools for analysis and plotting, there are great ways of integrating your code and plots into reports. You can do it the standard LaTeX way in `R` using Sweave() and Stangle (both part of the base distribution) to integrate with LaTeX (and make PDFs etc). If you want to use other formats (such as the widely used markdown format I use) then the `knitr` library is essential.

## Writing your own R packages
Roxygen, Shiny,...

If you are using R-studio, check out
https://support.rstudio.com/hc/en-us/sections/200107586-Using-RStudio

for many helpful suggestions.


## Getting more advanced with `R`?
Not surprisingly, Hadley Wickham (and others) have written several important books, including:
  - [Advanced R](http://adv-r.had.co.nz/) focused on really understanding how R works and best programming practices.
  - [R packages](http://r-pkgs.had.co.nz/) on how to write and build `R` packages.
  
  - Some may get a kick out of John Chambers book "Software for Data Analysis: Programming with R" on the evolution of `S` and `R`, and how to program in it.
  - [Rcpp](http://www.rcpp.org/) for a great way of integrating c++ code and calling it from `R`
  - [devtools]() is helpful for developing R libraries
  - [testthat](https://github.com/hadley/testthat) for units tests
  - [Efficient R Programming](https://csgillespie.github.io/efficientR/) by Colin Gillespie and Robin Lovelace. Looks like a great book.
