# Bio720 In class R exercise


# how do we navigate to the correct directory?

setwd("~/Dropbox/Bio720/Sept28_data/")

# import data into R
file1 <- read.table("M200_sm_male_hdhorn_ACAGTG_L008_results.xprs", h=T)
file2 <- read.table("M172_sm_male_hdhorn_ATCACG_L003_results.xprs", h=T)

file1[1:3,1:3]
file2[1:3,1:3]
str(file1)

file1_ordered <- file1[order(file1$target_id), ]
file2_ordered <- file2[order(file2$target_id), ]

file1_ordered[1:3,1:3]
file2_ordered[1:3,1:3]


# with(file1,
    # file1_ordered <- file1[order(target_id), ])
