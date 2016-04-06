# Introduction to `R` Activities for Sept. 28th.

Today in class you will work in groups of 2-3 to import some data, examine it, and do some simple data manipulations.
We will be working with some counts of sequence reads from an RNAseq experiment for a Rhinocerous beetle.

![alt text](https://upload.wikimedia.org/wikipedia/commons/thumb/a/a7/Kabutomushi-20070710.jpg/220px-Kabutomushi-20070710.jpg "Rhino beetle!!!")


The counts are from a program called [eXpress](http://bio.math.berkeley.edu/eXpress/overview.html). As we will examine in class, this generates a number of different count statistics. We will talk more about them later in the class. However, we need to decide which measures of transcript abundance to use (for now we will just use TPM).

Among the issues we have to deal with, the "target_id" names represent the current computationally generated (partially based on blast matches to another closely related beetle) names for genes. We may want to fix those. More importantly, `eXpress` does not sort the genes in any sensible way. So the order of genes in one file (from one sample) is different than the order in another file.

So we will want to:

1. Download the files
2. Sort the files by the identifier (target_id) so each file is in the same order.
3. Extract each of the columns for the TPM measure for the 4 files, and place them in a single data frame or matrix.
4. Look at correlations between measures of TPM between each of the samples.
5. Perform scatterplots for these expression measures.
6. Generate a subset of the "genes" which have an average TPM (across all samples) greater than 10.
7. Time permitting make an [MA plot](https://en.wikipedia.org/wiki/MA_plot). The two treatment groups are the "lg" (for large) VS. "sm" (for small).

## Download the files

I have made copies of the files you will need at
```bash
/2/scratch/bio720_ID
```
On the Golding lab server in which you each have accounts. You should navigate there and copy the files to your local machine.


Alternatively I have put them in a dropbox folder at
https://www.dropbox.com/sh/iko5nbtlg045b0q/AADPc5Fx2bQdexabn-2obsfua?dl=0
