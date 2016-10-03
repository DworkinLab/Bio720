
Here I am just providing some links to various `R` resources that may be helpful as you start to learn to program and perform analyses in `R`.

## R Tutorials
[R studio tutorials](https://www.rstudio.com/online-learning/)

[swirl interactive R tutorials](http://swirlstats.com/)

### Some of my favourite (general) R books
(less about understanding and performing analyses, more about programming and using R). There are many many more, but these are just a few I have used and read.

[R programming for bioinformatics](http://www.crcnetbase.com/isbn/9781420063684). Robert Gentleman. A very nice introduction to programming in R. Not much to do about bioinformatics *per se*, but a good way of learning `R`.

[The R Book](http://www.mcmu.eblib.com.libaccess.lib.mcmaster.ca/patron/FullRecord.aspx?p=297479). A really nice introduction to using R, data structures in R, plotting and basic statistics. If you are going to be using R, you want to have a copy of this to use!

[The art of R programming](https://www.nostarch.com/artofr.htm). A really nice introduction to programming in R. Does not assume you have any background in programming.

[R in a nutshell](http://shop.oreilly.com/product/0636920022008.do). Once you have some basic understanding of programming in general, this is a nice place to look things up quickly and does a good job of explaining what R is doing with respect to data structures.


## how to find libraries that you might need?
[CHECK OUT CRAN task views](https://cran.r-project.org/web/views/) to see links to libraries with specialized functions.

For genomics and bioinformatics most R libraries are available through [bioconductor](http://bioconductor.org/). There are hundreds of amazing libraries. One thing to know is that some of the way that data is stored for some bioconductor objects may seem a bit different than other R objects. It takes a bit to get used to it, but they can be very helpful down the road!

*The Hadley Wickham-verse*. [Hadley Wickham](http://hadley.nz/) has developed some very important and widely used libraries in R, and is one of the chief scientific officers for R-studio. In particular one of his goals was to make doing routine things (data munging, aggregating, cleaning and plotting) in R simpler and more intuitive than is found with standard R syntax, and as such many of his libraries are widely used (https://github.com/hadley). It is worth noting that of his effort to do so, he has created a sort of macro language in R (more like a parallel syntax really) that can be quite different. However, while it can be quite different than base R syntax, the learning curve is less steep (part of the point). These packages include: dplyr, reshape2, tidyr, readr, devtools, stringr, ggplot2, lubridate. All of these are available from his [github](https://github.com/hadley/) page. He has also written several books that cover the gamut of library specific things (ggplot2) to advanced R programming and writing your own packages.

data.table has some amazingly useful functions to import and work with large data.

Also check out [this list](https://support.rstudio.com/hc/en-us/articles/201057987-Quick-list-of-useful-R-packages) of other suggested packages (some are the same).

Making reports
In addition to great tools for analysis and plotting, there are great ways of integrating your code and plots into reports. You can do it the standard LaTeX way in `R` using Sweave() and Stangle (both part of the base distribution) to integrate with LaTeX (and make PDFs etc). If you want to use other formats (such as the widely used markdown format I use) then the `knitr` library is essential.

Writing your own R packages
Roxygen, Shiny,...

If you are using R-studio, check out
https://support.rstudio.com/hc/en-us/sections/200107586-Using-RStudio

for many helpful suggestions.


## Getting more advanced with `R`?
Not surprisingly, Hadley Wickham has written several important books, including:
  - [Advanced R](http://adv-r.had.co.nz/) focused on really understanding how R works and best programming practices.
  - [R packages](http://r-pkgs.had.co.nz/) on how to write and build `R` packages.
  - Some may get a kick out of John Chambers book "Software for Data Analysis: Programming with R" on the evolution of `S` and `R`, and how to program in it.
