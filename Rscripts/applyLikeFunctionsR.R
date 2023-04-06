# Using the apply functions - last updated Sept 22, 2015

# There are a set of really useful functions to repeatedly perform the same action (i.e. calling a function across the whole data set)

# apply, tapply,lapply, sapply (friendly version of lapply),mapply, rapply(recursive version of lapply) by, aggregate, sweep, replicate (wrapper for sapply)



# import and read in data


## YOU ONLY NEED TO DO ONE OR THE OTHER. IT IS THE SAME DATA AT TWO DIFFERENT LINKS

dll_data = read.csv("https://raw.githubusercontent.com/DworkinLab/DworkinLab.github.io/master/dataSets/Dworkin2005_ED/dll.csv", 
    header = TRUE, 
    stringsAsFactors = TRUE)

## OR ##

dll_data <- read.csv("http://beaconcourse.pbworks.com/w/file/fetch/35183279/dll.csv", 
    header = TRUE,
    stringsAsFactors = TRUE)   #data frame input


# Now we will change some of the variables to factors in R.
dll_data$temp <- factor(dll_data$temp)
dll_data$replicate <- factor(dll_data$replicate)

dll_data <- na.omit(dll_data) # normally I would not do this here, but this is just for demonstration, so I got rid of all rows with missing data



# let's take a look at the data
str(dll_data)

# let us look at the summary by column (MARGIN=2)
apply(dll_data, MARGIN = 2, summary) # MARGIN = 1 would be by row
# Not super useful, but you can see how we have some summary data of column
# Clearly this is not the most useful application of it since summary(dll_data) is more informative


# Using only the numeric columns
apply(dll_data[,5:8], MARGIN = 2, summary) # this gives us more useful information
apply(dll_data[,5:8], MARGIN= 2 , sd)  # standard deviations
#  apply is designed for matrices and arrays, and it coerces data frames ( as we used above) so be careful with data frames.
  

# what if we wanted to compute custom functions, such as the standard error of the mean for each column
# Since there is no built in function, we can write our own, and include it in the apply function

# Calculating the standard error of the mean
sem <- function(x){
	n <- length(x) # number of observations
	sd(x)/sqrt(n)  
	}
	
cv <- function(x) {
	sd(x)/mean(x)}

# Let's now apply these functions (using the apply() ) to a few columns of our data set.
	
apply(dll_data[,5:8], 2, sem)
apply(dll_data[,5:8], 2, cv)

# This approach can be used for any function you can think of.

# For people doing association mapping, or other genomic studies, this can be fruitfully utilized for going through the data. HOWEVER.. The apply family of functions dynamically grow the output vector/matrix, which while fast for small data sets, is terrible for really large data sets. In these situations it is best to pre-allocate memory for your output object, and fill it via a for loop.


# Often we are not interested in computing some function across rows or columns but across some treatment levels (i.e. to look at means across treatments for instance, compute a test statistic for each gene in the genome, or estimated slope for biomass~ temp for each lake in our study). 

# The "traditional" way is to use a for loop, which would look like this::

# Determine number of levels of line (line is an another term for a genetic strain)
N <- nlevels(dll_data$line)

# Generate an empty vector of length N
LinesSctMeans <- rep(NA, N)

#Set the names for each row of your empty vector with the names of the levels of the factor. This can be important for clarity and for running it (iterating across a vector of factor levels).
names(LinesSctMeans) <- levels(dll_data$line)

# Then the actual for loop:
for (i in levels(dll_data$line)){
	LinesSctMeans[i] <- mean(dll_data$SCT[dll_data$line==i])
}

LinesSctMeans

#####

#There is a much easier way to do this using a useful function called tapply, which applies a function along factor levels::

with(dll_data, tapply( X = SCT, INDEX = line, FUN = mean))



# We can also "cross" factors such as temp and genotype to get the means for all 4 "cells"
with( dll_data, tapply(SCT, temp:genotype, mean))

# You may instead prefer to have this in a 2x2 table format. This is easy to do. You just specify a list of factors. i.e. list(temp, genotype), instead of looking at their "interaction" (temp:genotype)

with( dll_data, tapply(SCT, list(temp, genotype), mean))

# Of course this is more useful for data sets for tables of more than 2x2


# As with apply, we can write our own functions and use them in tapply
cv.SCT <- with( dll_data, 
    tapply(SCT, list(line, genotype), cv))

cv.SCT

plot(cv.SCT, 
    xlim = c(0.05, 0.25), ylim = c(0.05,0.25), 
    xlab = " CV for Dll mutant allele",
    ylab = "CV for wild type Dll allele",
    pch = 16) 
abline(a = 0, b = 1)

# This result is actually pretty cool if you think about it!
	
	
# You should think about how you want the output represented. For example...	
( X.1 <- with( dll_data, tapply( SCT, list( line, genotype), cv)) )
( X.2 <- with( dll_data, tapply( SCT, line:genotype, cv)) )

dim(X.1) # 27 rows by 2 columns
dim(X.2) # one vector of length 54. 



# A really useful application of tapply. Computing Levene's statistic for multifactorial designs.

# One method to compare patterns of variance (i.e heterogeneity in the variances of different treatments) is to use Levene's Statistic. Indeed a simulation study by Schultz (1985) advocated the use of log transformed variables, and the median as a measure of central tendancy.  You can often find Levene's statistic functions for the equivalent of one-way ANOVA like designs, but it is less common for multi-factiorial designs that include crossed factors (see Dworkin 2005 Evolution).

# I originally wrote code for this is SAS, but I found a function in R for Levene's Statistic to demonstrate how easily it can be done (got it from the web, I wish I knew its original source). The version below is my modification of the original function for multi-factorial designs (sensu Dworkin 2005)

attach(dll_data)

#original version, designed for a single group (unknown source for function)
levene.test <- function(y, group) {
    meds <- tapply(y, group, median, na.rm = TRUE)
    resp <- abs(y - meds[group])
    table <- anova(lm(resp ~ group, na.action = na.omit))
    row.names(table)[2] <- " "
    cat("Levene's Test for Homogeneity of Variance\n\n")
    table[,c(1,4,5)]
    }

# work in progress - but it works for multiple groups 
levene.test <- function(y, group) {
    meds <- tapply(y, group, median, na.rm = TRUE)
    abs(y - meds[group])}
    
lev_stat <- with(dll_data, levene.test(SCT, genotype:temp:line))
ls.lm <- lm(lev_stat ~ genotype*temp*line, data=dll_data)
anova(ls.lm)  # of course this ANOVA is not particularly valid given the lack of normality, use bootstrapping or permutation.


########## General Purpose Levene's Deviates Calculator
### Below I have written a general purpose function to generate Levene's deviates
LeveneDeviates <- function(y, group, med=TRUE, log_trans=TRUE) {
    
    #log transform data?
    if (log_trans==TRUE)
        y = log(y)
    
    # allows for use of mean or median as measure of central tendency
    if (med==TRUE)
        meds <- tapply(y, group, median, na.rm=TRUE)
    else 
        meds <- tapply(y, group, mean, na.rm=TRUE) 
    
    # calculates deviates for each observation from a measure of central tendency for a "group"
    abs(y - meds[group])}
    
lev_stat <- with(dll_data, 
    LeveneDeviates(y = SCT, group = genotype:temp:line, med=TRUE, log_trans=FALSE))

ls.lm <- lm(lev_stat ~ genotype*temp*line, data=dll_data)
anova(ls.lm)  # of course this ANOVA is not particularly valid given the lack of normality, use bootstrapping or permutation.

########    
    
#################    
    # by()
# The by function is a wrapper for tapply that can be useful. 

# For instance  if you are trying to run the same linear model over and over again, but with different subsets of data (like in RNA seq analysis, by gene. BUT DON"T USE THIS FOR THAT. THis is not a statistically sensible way! Just for demonstration)

# For this example, I am running a model on the effect of genotype on each line in the data set
model.extract <- by(dll_data, INDICES = dll_data$line, function (x) lm (SCT ~ genotype, data = x))

# Looking at the object..
model.extract # We get the model info for each linear model
dim(model.extract) # 27 models in here

# If you just wanted to grab the slopes
sapply(model.extract, coef)[2,]

# Or extract other things
sapply(model.extract, logLik)  
    
    
#You can also do this for combinations of factors for models
model.extract.GxT <- by(dll_data, INDICES=list(dll_data$genotype, dll_data$temp), 
    function (x) lm (SCT ~ line, data=x))
    
dim(model.extract.GxT)
sapply(model.extract.GxT, coef)  




# It is important to note that while the apply functions are really convenient for small to moderate sized data tables ( <1000 factor levels or so), they can be slower than for loops for larger data sets (where you are applying across more than 1000 column/rows/factor levels). So in particular for many genomic studies, it is better to generate (pre-allocate) an empty matrix for your results and populate it with a for loop.


#####################    
# An example of the sapply function to use one of the pre-built power functions.

power = seq(from=0.1, to=0.9, by=0.01) # creates a vector of values from 0.1 to 0.9

# This function is just so that we can compute the sample size we need for varying levels of power
pow.test <- function(x){
	pow2 <- power.t.test(delta=0.5, sd=2, sig.level=0.05, power=x)
	return(pow2$n) # This pulls out the sample size we need
	}

power.n <- sapply( power, pow.test) # This just uses one of the apply functions to repeat the function (pow.test) for each element of the vector "power"

plot(power.n ~ power)
    
 #lapply will do the same as sapply, but outputs a list.
 
 

 
#replicate() see my examples in my resampling and power analysis scripts.  Very useful for Monte Carlo simulations or resampling where there is no dependencies on previous iterations.

# Say we wanted to simulate 100 random variables from a Normal Distribution with mean = 5, and sd = 2  (~N(5,2))
Sim.Normal <- replicate(n = 4, rnorm(100, 5, 2)) # 
