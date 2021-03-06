## ------------------------------------------------------------------------
x <- 1

for (i in 1:9) {
    x <- x + 1
    print(x)
}

x


## ------------------------------------------------------------------------
rm(x)

countFun <-function(x)
{
  
  for (i in 1:10) {
    x <- x + 1
    print(x)
}}

countFun(1)
x

## ------------------------------------------------------------------------
# initialize a vector of length 1
x <- NA

#system.time( YOUR FOR LOOP HERE)

system.time( for (i in 1:10000) {x[i] <- i^2})

## ------------------------------------------------------------------------
# initialize a vector of length 1
rm(x)
x <- rep(NA, 10000)

#system.time( YOUR FOR LOOP HERE)

system.time( for (i in 1:10000) {x[i] <- i^2})

## ------------------------------------------------------------------------
rm(x)
x <- 1:10000
system.time(x <- x^2)

## ------------------------------------------------------------------------
n <- 10000 # just run this once, each time you change n

rm(x)
system.time(x <- replicate(n = n, expr = rnorm(1, mean = 5, sd = 5)))


## ------------------------------------------------------------------------
rm(x)
x <- rep(NA, n)
system.time(for (i in 1:n) x[i] <- rnorm(1, mean = 5, sd = 5))

## ------------------------------------------------------------------------
rm(x)
system.time(x <- rnorm(n, mean =5 , sd = 5))

## ------------------------------------------------------------------------
rm(list=ls())

# sample size, n
n <- 20

# intercept, a
a = 3

# slope, b
b = 0.25

# independent/explanatory variable
x <- rnorm( n = n, mean = 10, sd = 2)

# response/dependent variable
y <- rnorm(n = length(x), mean = a + b*x, sd = 1 )

plot(y ~ x, pch = 20, cex = 1.5)

# use the `lm` function to fit a linear model (including  regression)
mod_1 <- lm(y ~ x)

summary(mod_1) # just to look at.

(p_val <- summary(mod_1)$coef[2,4])


## ------------------------------------------------------------------------

rm(list=ls())

nsims = 1000

# sample size, n
n <- 20

# intercept, a
a = 3

# slope, b
b = 0.25

pvals <- rep(NA, nsims)

for (p in 1:nsims) {
	x <- rnorm( n = n, mean = 10, sd = 2)
    y <- rnorm(n = length(x), mean = a + b*x, sd = 1 )
    mod_1 <- lm(y ~ x)
    pvals[p] <- summary(mod_1)$coef[2,4]}

# proportion less than 0.05?
sum(pvals < 0.05)/length(pvals) # estimate of the power.

hist(pvals)   

## ------------------------------------------------------------------------
rm(list=ls())
nsims = 1000

pow_sim_1 <- function(a, b, n,  std_dev){
	x <- rnorm( n = n, mean = 10, sd = 2)
    y <- rnorm(n = length(x), mean = a + b*x, sd = std_dev )
    mod_1 <- lm(y ~ x)
    return(summary(mod_1)$coef[2,4])
}

pvals <- replicate(n = nsims, pow_sim_1(a = 3, b = 0.25, n = 20, std_dev =1 ))

length(pvals)
sum(pvals < 0.05)/length(pvals)
hist(pvals)

## ------------------------------------------------------------------------
rm(list=ls())
N = 200  # Number of simulations for inner loop. You generally want this to be >1000. 

p = rep(NA, N) # initializing the vector to store the p values in the inner for loop. 

#Global Parameter values
a = 3 # intercept
b <- seq(from=0, to=0.5, by=0.05)

sample_size <- seq(from = 10 ,to = 100 ,by = 5)  # Incremently increasing sample size 

power.size <- numeric(length(sample_size)) # initializing the vector to store the power at each sample size for the outer for loop.

### initialize the matrix to store all of the power estimates
power.b <- matrix(NA,length(sample_size),length(b))

## Now the actual for loop
system.time(
for (k in 1:length(b))  # across the different values for the slope
 {
  
  b_b <- b[k]
  
   for (j in 1:length(sample_size))  # looping through the different sample_sizes

    {
   
      s_s = sample_size[j]
      for (i in 1:N)
      {
       x <- rnorm(s_s, mean = 10, sd = 2)  # simulate values of predictor
       y_det <- a + b_b*x             # deterministic part of model
       y_sim <- rnorm(s_s, mean = y_det, sd = 2)  # Simulate y|x values
       lm1 <- lm(y_sim ~ x)                    # fit model given simulation 
       p[i] <- coef(summary(lm1))[2,4] # You may want to extract a different p-value from the model.
	  
     }
    
      power.size[j] <- length(p[p<0.05])/N   # How many p-values are less than alpha (0.05)
   }
   
    power.b[,k] <- power.size 
}
)


# Now we graph it.
par(mfrow=c(1,2))

#3D perspective plot
persp(y=b, x=sample_size, z=power.b, col="blue", theta=-65, 
    shade=0.75, ltheta=45, ylab="slope", xlab="Sample Size", 
    lphi=30, zlim=c(0,1.25), ticktype = "detailed")

# contour plot
contour(z=power.b, x=sample_size, y=b, col="blue",  ylab="slope", xlab="Sample Size")

#fancy contour
filled.contour(z=power.b, x=sample_size, y=b, 
    ylim=c(min(b),max(b)), xlim=c(min(sample_size), max(sample_size)), 
    xlab="Sample Size", ylab="slope", color = topo.colors)

## ------------------------------------------------------------------------
expand.grid(letters[1:5], 10:6)

## ------------------------------------------------------------------------
#Global Parameter values
a = 3 # intercept
b <- seq(from = 0, to = 0.5, by = 0.05) # values for the slope

std_dev = 2
sample_size <- seq(from = 10, to = 100, by = 5)


# use expand grid to get all combinations of b (slope) and sample_size
b_N <- expand.grid(b, sample_size)  # may be worth looking at b_N to see what is being stored.
dim(b_N)
colnames(b_N) <- c("b", "sample_size") 

# Here is the function to generate the simulation and fit the model given the simulated data.
SimulatePower <- function(sample_size, b_b, a, std_dev){
	x <- rnorm(sample_size, mean = 10 , sd = 2)
	y_det <- a + b_b*x
    y_sim <- rnorm(sample_size, mean = y_det, sd = std_dev)
    lm1 <- lm(y_sim~x)
    pval <- coef(summary(lm1))[2,4]}

# We can use this for one sample_size and slope
check_it_works <- replicate(1000, SimulatePower(sample_size = 15, b_b = 0, a = 3, std_dev = 2))

hist(check_it_works, freq=T)


## ------------------------------------------------------------------------
p_values <- mapply(SimulatePower, 
    sample_size  = b_N$sample_size, b_b = b_N$b, 
    MoreArgs=list(a = 3, std_dev = 2)) 

## ------------------------------------------------------------------------
system.time(
rep_p <- replicate(n = 200, mapply(SimulatePower, 
    sample_size  = b_N$sample_size, b_b = b_N$b, # arguments to vectorize across
    MoreArgs=list(a=0, std_dev=1)) )   # other parameters
)

## ------------------------------------------------------------------------
dim(rep_p)

## ------------------------------------------------------------------------
power_lev <- apply(rep_p, MARGIN=1, 
    function(x) length(x[x <= 0.05])/length(x)) # how many p-values are less than 0.05

## ------------------------------------------------------------------------
grid_matrix <- matrix(data=power_lev, nrow=length(b), ncol=length(sample_size))

## ------------------------------------------------------------------------
persp(x=b ,y=sample_size, z=grid_matrix, col="blue", 
    theta=-10, shade=0.75, phi=15, d=2, r=0.1,
    xlab="slope", ylab="sample size", zlab="power", 
    ticktype="detailed")

filled.contour(z=grid_matrix, y=sample_size, x=b, 
    xlim=c(min(b),max(b)), ylim=c(min(sample_size), max(sample_size)), 
    ylab="Sample Size", xlab="slope", color = topo.colors)

