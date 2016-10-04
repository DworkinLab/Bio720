## ------------------------------------------------------------------------
x <- 1

## ------------------------------------------------------------------------
str(x)
mode(x)
typeof(x)

## ------------------------------------------------------------------------
y = c(3, 4, 5)

## ------------------------------------------------------------------------
z = 3:5

## ------------------------------------------------------------------------
mode(y)
mode(z)
typeof(y)
typeof(z)

## ------------------------------------------------------------------------
y == z
all.equal(y, z)
# BUT
identical(y, z)

## ------------------------------------------------------------------------
typeof(z)
# because of how division is computed in R, it converts to double
typeof(z/z)
 
#  though, addition, substraction and multiplication do not require this, so stay integer
typeof(z*z)

## ------------------------------------------------------------------------
rm( x, y, z) # clean up
x = 1
y = "1"
z = "one"
a = TRUE
b = "TRUE"

## ------------------------------------------------------------------------
x == y
all.equal(x, y)
identical(x, y)

## ------------------------------------------------------------------------
y == z
all.equal(y, z)
identical(y, z)

# Which might be obvious by
mode(y)
mode(z)

## ------------------------------------------------------------------------
y == z
all.equal(y, z)
identical(y, z)

# Which might be obvious by
mode(y)
mode(z)

## ------------------------------------------------------------------------
a == b
all.equal(a, b)
identical(a, b)

# Which might be obvious
mode(a); typeof(a)
mode(b); typeof(b)

## ------------------------------------------------------------------------
b1 <- as.logical(b)
is.logical(b1); is.character(b1)
b == b1
a == b
typeof(b1); mode(b1)

## ------------------------------------------------------------------------
a1 <- as.character(a)
is.character(a1)
a == a1
a == b
typeof(a1); mode(a1)

## ------------------------------------------------------------------------
y1 <- as.numeric(y)
y == y1
identical(y, y1)

x == y1
identical(x, y1)
typeof(y1); mode(y1)

## ------------------------------------------------------------------------
rm(x, y, z, a, b)

## ------------------------------------------------------------------------
x = T
x
y = TRUE
y

x == y
identical(x, y)

## ------------------------------------------------------------------------
sum(x)
as.numeric(x)

a = F
sum(a)

## ------------------------------------------------------------------------
sum(c(rep(T, 10), rep(F, 18), rep(T, 10), rep(F, 6)))

## ------------------------------------------------------------------------
ls()

## ------------------------------------------------------------------------
rm(list=ls())
ls()

## ------------------------------------------------------------------------
gene1 <- c(3, 4, 7, 9, 12, 6)
gene2 <- c(11, 17, 12, 25, 23, 7)
gene3 <- c(100, 103, 97, 94, 106, 111)

## ------------------------------------------------------------------------
genotype <- c("wildtype", "wildtype", "wildtype", "mutant", "mutant", "mutant")
genotype
mode(genotype)

## ------------------------------------------------------------------------
genotype2 <- rep(c("wildtype", "mutant"), each = 3)
genotype2
mode(genotype2)
class(genotype2)

## ------------------------------------------------------------------------
genotype3 <- gl(n = 2, k = 3, labels = c("wildtype", "mutant"))
genotype3
mode(genotype3)
class(genotype3)

## ------------------------------------------------------------------------
genotype2 == genotype3
identical(genotype2, genotype3)
all.equal(genotype2, genotype3)

## ------------------------------------------------------------------------
genotype2_factor <- as.factor(genotype2)
class(genotype2_factor)
mode(genotype2_factor)
identical(genotype3, genotype2_factor)
genotype3 == genotype2_factor

## ------------------------------------------------------------------------
genotype3_character <- as.character(genotype3)
genotype3_character 
class(genotype3_character)
mode(genotype3_character)
identical(genotype3_character, genotype2)
genotype3_character == genotype2

## ------------------------------------------------------------------------
day <- gl(n = 2, k = 1 , length = 6, labels = c(3, 6))
day
class(day)
mode(day)
typeof(day)

## ------------------------------------------------------------------------
as.character(day)

## ------------------------------------------------------------------------
as.numeric(day)

## ------------------------------------------------------------------------
as.numeric(as.character(day))

## ------------------------------------------------------------------------
gene_mat1 <- cbind(gene1, gene2, gene3)
gene_mat1

gene_mat2 <- matrix(c(gene1, gene2, gene3), nrow =6, ncol =3, byrow =FALSE)

## ------------------------------------------------------------------------
gene_mat1 == gene_mat2
identical(gene_mat1, gene_mat2)

## ------------------------------------------------------------------------
colnames(gene_mat2) <- c("gene1", "gene2", "gene3")
gene_mat2
identical(gene_mat1, gene_mat2)

## ------------------------------------------------------------------------
genotype <- rep(c("wildtype", "mutant"), each = 3)
day <- rep(c("3", "6"), times = 3)

genotype
day

treatments <- cbind(genotype, day)
class(treatments)
mode(treatments)

## ------------------------------------------------------------------------

all_the_data <- cbind(gene1, gene2, gene3, genotype, day)
class(all_the_data)
mode(all_the_data)

## ------------------------------------------------------------------------
rm(all_the_data)
all_the_data <- data.frame(gene1, gene2, gene3, genotype, day)
str(all_the_data)
class(all_the_data)
mode(all_the_data)


## ------------------------------------------------------------------------
all_the_data[ ,c(2:4)]


## ------------------------------------------------------------------------
all_the_data[c("gene2", "gene3", "genotype")]


## ------------------------------------------------------------------------
all_the_data$gene2; all_the_data$gene3; all_the_data$genotype


## ------------------------------------------------------------------------
all_the_data[["gene2"]]
all_the_data$gene2

## ------------------------------------------------------------------------
subset(all_the_data, select = c("gene2", "gene3", "genotype"))


## ------------------------------------------------------------------------
all_the_data$gene4 <- c(10, 11, 7, 11, 2, 3)

all_the_data

str(all_the_data)

## ------------------------------------------------------------------------
list_the_data = list(gene1, gene2, gene3, genotype, day)
list_the_data

str(list_the_data)
names(list_the_data)


## ------------------------------------------------------------------------
list_the_data = list(gene1 = gene1, gene2 = gene2, gene3 = gene3, genotype = genotype, day = day)
list_the_data

str(list_the_data)
names(list_the_data)


## ------------------------------------------------------------------------
list_the_data$random_variable = c(T,T,F) 

list_the_data

str(list_the_data)

## ------------------------------------------------------------------------
list_the_data$gene1

list_the_data[1]
list_the_data["gene1"]

list_the_data[[1]]
list_the_data[["gene1"]]

## ------------------------------------------------------------------------
class(list_the_data$gene1)
class(list_the_data[1])

class(list_the_data["gene1"])
class(list_the_data[[1]])

class(list_the_data[["gene1"]])

str(list_the_data$gene1)

str(list_the_data[1])
str(list_the_data["gene1"])

str(list_the_data[[1]])
str(list_the_data[["gene1"]])

## ------------------------------------------------------------------------
as.numeric(list_the_data[1])

## ------------------------------------------------------------------------
str(as.numeric(unlist(list_the_data[1])))

