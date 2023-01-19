###### Plotting a big PCA plot with scatter plot
# using RNAseqQC

library("RNAseqQC")
library("DESeq2")
library("ensembldb")
library("dplyr")
library("ggplot2")
library("purrr")
library("tidyr")
library("tibble")
library("magrittr")


# Load in data
sample_info <- read.csv("./plot_sample_info.csv")

# remove the outlier to match the counts data
sample_info <- sample_info[-5,]

salmon_counts <- read.csv("./salmon_counts_raw_all.csv", row.names = 1)

salmon_counts <- round(salmon_counts)

salmon_counts_mat <- as.matrix(salmon_counts)

# make it so that they match explicitly
rownames(sample_info) <- colnames(salmon_counts)

# make DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = salmon_counts, 
                              colData = sample_info, 
                              design = ~ background + perturbation)

# filter low count genes 
dds <- filter_genes(dds, min_count = 5)

# variance stabilizing
vsd <- vst(salmon_counts_mat)
# check if it worked! Think so
mean_sd_plot(vsd)

options(ggplot2.continuous.colour = "viridis")

# plot individual PC plots 
plot_pca(vsd, PC_x = 1, PC_y = 2, 
         color_by = "perturbation", 
         shape_by = "background", 
         point_alpha = 1, point_rel_size = 4)


# Plot scatter plot comparing PC1 through PC4
plot_pca_scatters(vsd, n_PCs = 4, 
                  color_by = "perturbation", 
                  shape_by = "background", 
                  point_alpha = 1, show_var_exp = FALSE)



