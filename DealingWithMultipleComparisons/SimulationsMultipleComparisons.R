# A little simulation to demonstrate why we need to correct for multiple comparisons.


#Let's assume we want to use an alpha of 0.05. That is 1/20 we are willing to reject P(H0|D).

# First we look at the expression of one gene. As it turns out gene1 is not differently expressed between our two treatments (males and females).  A better way of thinking about it, is that the expression of gene1 comes from the same underlying distribution in males and females.

# we go ahead and measure gene expression on 10 males and females.

males <- rnorm(10)
females <- rnorm(10)

dat <- data.frame(gene1 = c(males, females), 
                  sex = gl(2, 5, labels = c("M", "F")))
#we go ahead and run a t-test

t.test(gene1~sex, data=dat)

# We go ahead and now look at 1000 genes. As it turns out, in this tissue we are studying NO genes are differentially expressed (only we don't know that)


test.gene <- function(){
males <- rnorm(10)
females <- rnorm(10)

dat <- data.frame(gene1 = c(males, females), 
                  sex = gl(2, 5, labels = c("M", "F")))
#we go ahead and run a t-test

pval <- t.test(gene1~sex, data=dat)$p.value
return(pval)
}

pvals_1000genes <- replicate(1000, test.gene())

length(pvals_1000genes[pvals_1000genes < 0.05])

hist(pvals_1000genes)

# So what happened?


hist(runif(1000))  
# Just by chance you expect to get alpha% (that's the point). 

# So what do we do when we are looking at 1000 genes

# Can we adjust the p-values?

# There are many approaches to deal with this. The basic ones just adjust the p-value based on number of comparisons. So to "maintain" your nominal alpha=0.05 when you are doing 1000 comparisons, divide your alpha by 1000  I.e. 0.05/1000

adjusted_alpha <- 0.05/1000
# Which is really small
length(pvals_1000genes[pvals_1000genes < adjusted_alpha]) # No significant hits

# The problem is, that such an approach is WAY TOO conservative, and you are not maintaining your nominal alpha at all!

# Another approach is the sequential Bonferroni, which adjusts for each comparison alpha/1000, alpha/999, alpha//998 etc... However, for large numbers of comparisons this is almost as overly conservative as Bonferroni

# So there are a variety of methods that actually adjust in somewhat more sensible ways. THis includes an important set of approaches based on False Discovery Rates (FDR) as opposed to false positives. Here we are not asking about how many false positives in the full set of comparisons, but how many "false positives" given a set of true positives. The most commonly used approaches in genomics are Benjamini Hochberg, as well as q-values. To a first approximation these methods use the differences between the theoretical expectation of the distribution of p-values (assuming nothing interesting) VS. observed to determine the false discovery rates (and adjustments)

# In R most of these approaches can be
p.adjust()

# For the q-value approach there is a seperate library in R.