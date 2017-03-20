## How about mixed models

Since line should be fit as a random effect, we need to consider our options. In R the best choice is to then go ahead and use `MCMCglmm` which has a natural interface for general multivariate mixed models.  It takes a while to get used to the interface, but here is an example.


```{r}
fmla.T2  <- as.formula(paste("cbind(femur_s, tibia_s, tarsus_s, SCT_s)" ,"~", "trait - 1 + genotype + temp + genotype:temp"))
fam.test <- rep("gaussian", 4 )


# Prior for model. 
prior.model.1 <- list( R=list(V=diag(4)/4, nu=0.004),  
                       G=list(G1=list(V=diag(4)/4, nu=0.004)))


# MCMCglmm
MMLM1.fit <- MCMCglmm(fmla.T2,
                      random=~ us(trait):line, 
                      rcov=~ us(trait):units,
                      prior=  prior.model.1,
                      data= dll_data, 
                      family = fam.test, nitt= 10000, burnin= 2000, thin=10)
```

Let's take a look.

```{r}
summary(MMLM1.fit)
```

