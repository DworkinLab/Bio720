---
title: "Bio720 - Simulating Data"
author: "Ian Dworkin"
date: "11/19/2018"
output:
  ioslides_presentation:
    incremental: yes
    keep_md: yes
    widescreen: yes
  slidy_presentation:
    incremental: yes
---



## Slide 1

One of the most important skills in your bag as new computational biologists is the ability to perform simulations. In particular, to simulate data or evaluate models numerically (or both).

## Slide 2
In biology mathematical models are the basis for much of the theoretical and conceptual background in disciplines like ecology, evolution, population genetics, molecular evolution & biochemistry.

## Slide 3

Many of the models that are developed are not *analytically* tractable. That is, without making very strong biological assumptions there are no *closed form solutions*.

That is the models can not be solved for the general case.

## Using computers to find numerical solutions

- For such models, "solutions" can still be found.
- For some models, stability analysis can be performed (among other techniques).
- However, most often, scientists resort to computers to identify *numerical solutions*

## Deterministic VS. Stochastic
- Within *dynamical* models, that include the majority of models in biology there are two broad categories.
- **Deterministic** - where the outcome of the model entirely depends on the model itself and the starting conditions.
- **Stochastic** - where random events influence the outcomes of the model, and only the probability of events can be predicted, not exact outcomes.

## A simple deterministic model one-locus model of natural selection.
- In population genetics we often start with a simple deterministic model of how selection on an allele influences its frequency in the population over time.
- In a simple haploid model we can envision two alternative alleles, *A* and *a*.
- *a* is initially fixed in the population, but the new mutation *A* arrives and it is beneficial. What will happen to this allele?

## Haploid selection model
- Frequency of *A* at time *t* is $p(t)$
- Fitness of *A* is $W_A$, and for *a* is $W_a$

$p(t+1) = \frac{p(t) W_A}{p(t) W_A + W_a(1-p(t))}$

$  = \frac{p(t) W_A}{\bar W }$

- Where $\bar W$ is mean fitness for the population

- How would you convert this into *R* code for one generation to the next?

- Start with what variables you need.
## Converting this into *R* - What variables
- We need variables for the fitness values for each allele, $W_A$, $W_a$ and for allele frequency of A $p(t+1)$ at time t and t+1.

- We can put the pieces together. Write an expression for mean population fitness called `w_bar`

$ \bar{W} = p(t) W_A + W_a(1-p(t))$




## Slide with Bullets

- Bullet 1
- Bullet 2
- Bullet 3

## Slide with R Output


```r
summary(cars)
```

```
##      speed           dist       
##  Min.   : 4.0   Min.   :  2.00  
##  1st Qu.:12.0   1st Qu.: 26.00  
##  Median :15.0   Median : 36.00  
##  Mean   :15.4   Mean   : 42.98  
##  3rd Qu.:19.0   3rd Qu.: 56.00  
##  Max.   :25.0   Max.   :120.00
```

## Slide with Plot

![](Bio720_SimulatingData_files/figure-html/pressure-1.png)<!-- -->

