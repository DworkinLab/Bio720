##### Last modified June 13th, 2011. Ian Dworkin
# Testing conditions.

 # The standard if else
p.test <- function(p){
    if (p <= 0.05) print ("yeah!!!!") 
        else if (p >= 0.9) print ("high!!!!")
        else print ("somewhere in the middle")
}

# Also the ifelse(test, yes, no) function. ifelse() is far more useful as it is vectorized.
p.test.2 <- function(p){
    ifelse(p <= 0.05, print("yippee"), print("bummer, man")) }


# Test this with the following sequence. See what happens if you use if vs. ifelse(). 
x <- runif(10, 0, 1)
 
 
# There are many times that you may think you need to use an if with (iterating with a for loop... see below), or ifelse, but there may be far better ways.

# For instance, say you are doing some simulations for a power analysis, and you want to know how often your simulation gives you a p-value less than 0.05.

p.1000 <- runif(n=1000, min = 0, max =1 ) # generate 1000 values between 0-1, which we will pretend are our p-values from our simulation.

# you may try and count how often it less than 0.05
p.ifelse <- ifelse(p.1000 < 0.05, 1, 0) # If it is less than 0.05, then you get a 1, otherwise 0. 

sum(p.ifelse)/length(p.1000) # Our approximate false positives. Should be close to 0.05

# However the best and fastest way to accomplish this is to use the index, by setting up the Boolean (TRUE/FALSE) in the index of the vector.
length(p.1000[p.1000 < 0.05])/length(p.1000) # Same number, faster and simpler computation. 



# Simple loops

# while() function..  I tend to avoid these, so you will not see them much in class......
i <- 1
while(i <= 10){
    print(i)
    i <- i + 0.5}

	
	
	
# for loop, automatically iterates across a list (in this case the sequence from 1:10)
for (i in 1:10){
	 print(i)}
	 

# if you do not want to use integers, how might you do it using the for() ?

for (i in seq(from=1,to=5, by=0.5)){
	 print(i)}


#### behavior of strings.	
# Using strings is a bit more involved in R, compared to other languages. For instance the following does not do what you want::

for (letter in "word") {
	print(letter)
}	
# (try letters for a hoot.)	

# Instead in R, we have to split the word "word" into single characters using strsplit(), i.e::

strsplit("word", split="")

# So for the for loop we would do the following::
for (letter in strsplit("word", split="")) {
	print(letter)
}	



# Traditional way of outputting random numbers in many programming languages

for (i in 1:100){
	print (rnorm(n=1, mean=0, sd=1))} # We are cycling through and generating one random number at each iteration.	 Look at the indices, and you can see we keep generating vectors of length 1. 


##############
# What if we wanted to put all of these numbers in a vector?

# First we initialize a vector to store all of the numbers. Why do we initialize this vector first?	
n <- 100000
x <- rep(NA,n) 
# The step above creates a vector of n NA's. They will be replaced sequentially with the random numbers as we generate them. 
head(x)

# equivalently R lets you express this like
# x <- numeric(length=n) # creating a vector with "n" 0's

# Now we run the for loop.
for (i in 1:n){
	x[i] <- rnorm(n=1, mean=0, sd=1) # for each i, one number is generated, and placed in x
	}	

head(x)	 

#####

# However this is computationally inefficient in R. Which has vectorized operations.

system.time(

for (i in 1:n){
	x[i] <- rnorm(n=1, mean=0, sd=1) # for each i, one number is generated, and placed in x
	})


# We can also use the replicate function to do the same thing. Easier syntax to write.
system.time(
	z <- replicate(n, rnorm(n=1, mean=0, sd=1))
	)
# ~ 20% faster on my computer	


# However, since R is vectorized, both of these will be far slower than	....
system.time(	
y <- rnorm(n,0,1)	
)

# About 65 times faster than the for loop

# The general rule in R is that loops are slower than the apply family of functions (for small to medium data sets)  which are slower than vectorized computations.