# Some behaviours of for loops


# looping through a vector
# First we have some data. Here are some Illumina HiSeq reads for one of our recent projects
read_1 <- "CGCGCAGTAGGGCACATGCCAGGTGTCCGCCACTTGGTGGGCACACAGCCGATGACGAACGGGCTCCTTGACTATAATCTGACCCGTTTGCGTTTGGGTGACCAGGGAGAACTGGTGCTCCTGC"

read_2 <- "AAAAAGCCAACCGAGAAATCCGCCAAGCCTGGCGACAAGAAGCCAGAGCAGAAGAAGACTGCTGCGGCTCCCGCTGCCGGCAAGAAGGAGGCTGCTCCCTCGGCTGCCAAGCCAGCTGCCGCTG"

read_3  <- "CAGCACGGACTGGGGCTTCTTGCCGGCGAGGACCTTCTTCTTGGCATCCTTGCTCTTGGCCTTGGCGGCCGCGGTCGTCTTTACGGCCGCGGGCTTCTTGGCAGCAGCACCGGCGGTCGCTGGC"

# q1, what species are they from?

reads <- c(read_1, read_2, read_3)


# say we wanted to print each character from read_1, how do we do this using a for loop

for (i in read_1) {
	print(i)
}

# While this prints out the full string, it does not do so character by character. 


# We can make this one string into a vector of strings with each element being a single letter

read_1_split <- strsplit(read_1, split = "", fixed = T) # a list

# Since it is a list, which is annoying we want to coerce this into a character vector so we can work with it more easily.

read_1_char <- as.character(unlist(read_1_split))
mode(read_1_char)

# Printing the list once.
for (i in read_1_split){
	print(i)
}

# VS looping through each character

for (i in read_1_char){
	print(i)
}


# How about if we wanted the number of occurrences of each base? Or their frequencies? Well, to start with the easiest thing to do would be to use the table() function!

# This works for lists
table(read_1_split)/lengths(read_1_split)  # note the use of lengths() not length()

# for a vector of characters you can use this
(table(read_1_char)/length(read_1_char))


# how would you make this into a nice looking function that can work on either lists or vectors of characters?


BaseFrequencies <- function(x) {
    if (mode(x) == "list") {
    	tab <- table(x)/lengths(x)}
    
    else {
    	tab <- table(x)/length(x)
    }	
   print(tab, digits = 2)
   print(mode(x), digits = 2)
}


# Add in what to do if it could also come in as a single multi-character string.


# Now try doing it with sapply across the three reads

# might be worth re-writing this more generally
BaseFrequencies <- function(x) {
    
    # if it is a single string still
    if (length(x) == 1 & mode(x) == "character") {
    	x <- strsplit(x, split = "", fixed = T) 
        x <- as.character(unlist(x))
    }     
    
    if (mode(x) == "list") {
    	tab <- table(x)/lengths(x)}
    
    else {
    	tab <- table(x)/length(x)
    }	
   print(tab, digits = 2)
   print(mode(x), digits = 2)
}

sapply( reads, BaseFrequencies, USE.NAMES = F)
