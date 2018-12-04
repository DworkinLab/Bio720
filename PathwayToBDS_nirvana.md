# Pathway to Bioinformatics Data Skill Nirvana

## Introduction
You have done your readings through sections of Vince Buffalo's Bioinformatics Data Skills book. You have done DataCamp tutorials and developed some command line skills in UNIX and the basics of a programming language like R or Python. You have seen various command line tools (like bwa or samtools), and many R libraries. You have done all the class activities. Now it is time for you to do your own research....

 The only question is how do you begin? There are seemingly so many command line tools to use, or programming languages that you might want to use like R, Python or Perl. Within each of those languages you need to choose whether to code everything yourself with base functions and DIY, or use established libraries (like *stringr* or *stringi*), and even very specific libraries for dealing with genomics data (*Biostrings* in bioconductor/R). It is all very overwhelming and difficult to get started.

As such here are a few thoughts. Note, these thoughts are just my own opinions based upon my own learning, watching my students learn.  Also, while I give specific examples of tools, these are not meant as recommendations for "best". Indeed I am largely arguing that most tools that you find after a bit of google-fu and reading will serve you, and it pays to spend time learning a few tools well initially for your pipeline. You don't need to know every command line program or *R* library out there!


## Be patient with yourself and the software tools, especially at the beginning.

Hopefully as your skills have developed in class, you have seen that things that took you several hours back in September, only take a few minutes now. That is fantastic. Keep this in mind as you move forward. Getting started with new analyses, new software or a new programming language with take time, with challenges in knowledge and skill development along the way. That is perfectly normal, and the best advice is patience, resilience and practice as you learn. If you have a supervisor (i.e. a PI) who does not regularly do genomic or bioinformatic analyses you may want to send them this link, as they might think this should be as easy as using an Excel spreadsheet.

The good news is (and you probably have already got a sense of this), as you develop some basic UNIX scripting skills and pick up a programming language, learning how new languages, libraries or tools work becomes easier and more familiar, including the language of help or man pages. So [keep on keepin' on!](https://youtu.be/4QDECwqKE0g?t=13)


## choose a path

The fact is that there are many different ways to achieve the same goals for many common bioinformatic analyses. For instance, imagine you are identifying and counting sequence motifs in a set of genomes. You can absolutely accomplish this with command lines tools alone, or using Bioconductor/R, or Biopython. While there might be some very good reasons to pick one (command line tools that work on streaming data may make sense for very large genomes for instance). You will likely be able to accomplish 97.5% of the same tasks with any of these approaches. Pick one approach and learn it well to start with.

If you are going the command line route, spend your energy learning UNIX tools well, including programs like sed and awk for more powerful pattern matching/text processing and picking key, well used general bioinformatics tools that make sense for your work (BLAST, samtools, [seqtk](https://github.com/lh3/seqtk), [emboss](http://emboss.sourceforge.net/) and then use more specialized tools ([meme](http://meme-suite.org/index.html), [rast](http://tagc.univ-mrs.fr/rsa-tools/), [bedtools](https://bedtools.readthedocs.io/en/latest/)). It is usually easier to find tutorials and good working examples that you can base your own analysis off of from some of these more widely used tools. You will still need to learn the basics of a programming language like *R* or *python* for the data analysis and plotting, but you may only need them for those uses to start with. Your *R* skills (for instance) are probably already close to good enough to get you going from what we have done in class!

If you are going the R/Bioconductor route, get to know the basics first (Biostrings, iRanges, GenomicRanges, GenomicFeatures). Spend the time going through the tutorials. Bioconductor makes this very easy (videos on youtube, written tutorials, lots of examples). You would still need to learn the basic UNIX skills for mapping reads etc... But probably not much more than you have already learned in Bio720!

While I am less familiar with the offerings of bio-perl and Biopython, they absolutely have the same broad set of general tools that one could choose.

### Utilize your local expertise (your lab, department or University) in choosing your path.

In addition to picking a path that best suites your preferences for language/environment, also consider other "experts" nearby. What do other individuals in your lab, or nearby labs use? If there is some local expertise, utilize this to your advantage. This may shape your initial approach, but may also help you gain deeper insights as well. Then once you are more comfortable you can delve further and deeper.

To re-iterate all of this, don't try to learn all of the tools at once. Pick an approach (UNIX shell, R, python or perl), go deep with that, and for the moment only learn enough of the other skills to enable you to finish the analyses.

## DIY isn't always best. Don't re-invent the wheel!

One bad habit you could potentially pick up from the activities we have done in Bio720 is to write your own function to do *X*, or pipeline to do *Y*, even when there are many packages with functions that allow you to do this with a single line of code! Exercises in class are  motivated by pedagogical goals to reinforce concepts and consolidate knowledge. We rarely design them as a reasonable, practical (nor efficient) method to go about doing a particular task. You have already seen this with the DIY data munging activities, as opposed to using base R functions like colMeans, or even dplyr tools which can make it even easier.

The fact is that many of the libraries that are well used for common tasks have also (probably) been well vetted in their functionality. I could write some function with loops with `gregexep()` to  pull out and count some sequence motif across many sequences. Or, I could use the functionality from *stringr* or *stringi* libraries and accomplish the same task far more quickly, and be much less likely to introduce bugs into my code. Better yet, if I was working with nucleotide sequence data (where I need to think about both strands of DNA, including the reverse complement), then *Biostrings* may be the choice to start with. Indeed, *Biostrings* has its own functions (including its own version of `gregexpr`)that work when multiple instances of the pattern overlap with one another. No need for a DIY version! Ditto for command line tools. I could do it all with powerful (and general) tools like *sed* or *awk*, but maybe the tools in *meme* or *emboss* will do it more quickly, and they have been well vetted for bugs (and maybe even how fast they work) and nucleotide sequence data.

There will be plenty of functions that you will need to write DIY, but for basic operations, chances are there are plenty of functions or command line tools to use. So while writing your own functions is essential for many tasks, don't re-invent the wheel.


## Spend the time going through tutorials and examples.
  Many of the libraries and software packages have tutorials. Spend the time going through these to get comfortable. Don't be afraid to experiment with them, and then spend some time looking at the documentation to get a better sense of the parameters of the functions or programs and what can (and can't) be accomplished by it.  You don't always need to rush in with your own data, as this will often lead to a mess as you try hard to get the program to spit out any kind of output, even when it might be nonsensical.

## Avoid the GUI version of tools.
  While in genomics there are many GUI versions of tools like BLAST, or Galaxy for NGS data, etc.. I still recommend using your own code. While GUI tools can usually spit out logs (for the underlying code), researchers often forget to do this, and in the end forget how they did the analysis, what parameters they set, and even what versions of the software they used. Keeping your scripts is an important starting point for good reproducible research.

## Don't blindly trust the computer output

Computers are powerful, but exact. You can easily get unintended output because of a misunderstanding of what the program is intended for or capable of. This is doubly true for bioinformatic/genomic software tools. After hours of grinding through errors, sometimes you are so happy to get any output as an "answer" that you presume it is working. This assumption may not be warranted, so check! That is, you are the biologist, and have a sense of likely and possible values in the results. Do your results make sense? Can you check a little bit by "hand" to confirm you can get the same results? If not can you simulate some data with known properties to check?

This advice also extends to the software. Don't blindly trust it either. These software packages are being developed by researchers like you. In addition to accidental software bugs, in any analysis pipeline decisions are made as are assumptions. So two different software packages ostensibly performing the same process likely will have different answers. This can actually be an advantage. In our lab with some known biases of syntenic read mapping software (we used *bwa* in class as one tool for this), we only uses polymorphisms that we find in common when we run our pipeline identically across two different read mappers.

## Don't forget to visualize intermediate data
  In a rush to get the finish line of a complete analysis, researchers often skip some important QC & sanity checks along the way. Yes, most people do simple visualizations of reads in a `.fastq` file using `fastqc` before and after trimming. How about after you have mapped your reads and made your `.sam` or `.bam` files? Have you looked at your coverage across the genome? Not only some histograms, but actual coverage across sites?  How about after de-duplication. You don't expect this to change much, but have you checked to make sure? Ditto for removing poorly mapped reads, or repeat masking or when re-align indels in GATK?  It is important that for each step of your pipeline you:

  - write down a few points about what effect this should, what the output should be,.
  -  what you expect the output to look like, not just the file itself, but in terms of coverage genomically or locally (i.e. within a transcript for RNA-seq).
  - Possible pathologies you should check for.  

 While it may take you longer the first time with a new set of data (or using new tools), it will likely save you a lot more time of reaching the end, only to realize your analysis got horribly garbled along the way (and you have to re-trace your steps). Plus it may be very interesting in its own right. Maybe you will discover a big increase in local coverage suggesting a possible CNV that was not known, and you could check out! Remember to redo these checks if you introduce a new version of one of the programs or libraries into your pipeline!

I should also note (not to add more programs to learn) that there are often already great tools out there to vett these intermediate data files and look for known pathologies!

## keeping track of versions of software, libraries and genomes

  Since you want your research to be easily reproducible, knowing which versions of all the tools you used is important. So make sure that you add to your script some lines that automatically report this at the end of your analysis. In *R* this is easily accomplished. At the end of your script (`.R` or `.Rmd`) make sure to add a line like
```{R}
  sessionInfo()
```

Which will report the computer, version of R, and version of all R libraries.

If you have a UNIX shell script, most programs enable something like

```{bash}
  fastqc -version
```
Which you can list at the end of your script.


Ian Dworkin, with comments by Brian Golding and Ben Bolker
