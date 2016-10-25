set.seed(987654321)
# 1
StdErr <- function(x) {
	return(sd(x)/sqrt(length(x)))}

MyFunction <- function(x){
	return(c(mean = mean(x), 
	    median = median(x),
	    std_dev = sd(x),
	    stdErr = StdErr(x)))}

# 2
gene1 <- rnorm(mean = 10, sd =1, n = 1000)

MyFunction(gene1)

#3 
treatment1 <- gl(n = 2, k = 500, labels = c("treatment_A", "treatment_B"))

#4 
genotype <- gl(n = 2, k = 1, length = 1000, labels = c("wildtype", "mutant"))

#5
dose <- gl(n = 2, k = 20, length = 1000, labels = c("highdose", "lowdose"))

#6
gene2 <- rnorm(n = 1000, mean = 50, sd = 5)
gene3 <- rnorm(n = 1000, mean = 1000, sd = 50)

#7
gene_dat <- data.frame(treatment1, genotype, dose, gene1, gene2, gene3)

#8
MyFunction(gene_dat$gene1)  
MyFunction(gene_dat$gene2)
MyFunction(gene_dat$gene3)

# could also do:
MyFunction(gene_dat[,4])
MyFunction(gene_dat[["gene1"]])

## or more succicently
apply(gene_dat[,4:6], MARGIN = 2, MyFunction)

#9 - easiest way
MyFunction(gene_dat$gene1[gene_dat$gene1 < median(gene_dat$gene1)])
MyFunction(gene_dat$gene2[gene_dat$gene2 < median(gene_dat$gene2)]) 
MyFunction(gene_dat$gene3[gene_dat$gene3 < median(gene_dat$gene3)]) 

#10

# Cleanest to make a new helper function

SEM <- function(x) {
	up <- mean(x) + StdErr(x)
	down <- mean(x) - StdErr(x)
	return(c(up = up, down = down))
}

SEM1 <- SEM(gene_dat$gene1)
MyFunction(gene_dat$gene1[  gene_dat$gene1 < SEM1["down"] | gene_dat$gene1 > SEM1["up"]])
 
SEM2 <- SEM(gene_dat$gene2)
MyFunction(gene_dat$gene2[  gene_dat$gene2 < SEM2["down"] | gene_dat$gene2 > SEM2["up"]]) 

# etc