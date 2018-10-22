# Getting used to markdown


## What is markdown?

Markdown is a specific *markup* language (yes, a play on words) for generating documents (html and LaTeX are other examples). Markup languages  enable the generation of *dynamic documents* that allow for the generation of reports, presentations and manuscripts. Markdown (compared to html or LaTeX) is very simple and readable file format. that enables the generation of *dynamic documents* that allow for the generation of reports, presentations and even drafts of papers. Importantly it makes embedding blocks of code really easy, so it can help you if you want to keep an electronic copy of your lab notebook. As we will see soon (next week), with version control (like `git`) this can enable a powerful way of organizing your projects.  Since markdown is so simple, it really only takes about 5 minutes to get started, and you will be comfortable with it after a couple of days. There are many introductions to markdown on the web. I usually recommend [this one](https://daringfireball.net/projects/markdown/basics), although if you decide (after I introduce it) to use github, I suggest looking at [this as well](https://help.github.com/articles/markdown-basics/). Importantly for us, there is `R flavoured` markdown, which you can find more about [here](https://rmarkdown.rstudio.com/).

The goal of `markdown` is to keep everything human readable even as it allows for some simple formatting and generation of figures that can be imbedded. This differs from other *markup* languages While the markdown (using the extension `.md`) files are just regular text files, they can be converted to PDF, html and many other formats easily (including .tex). Indeed a number of text editors (google "markdown editor") can render the markdown in realtime. Otherwise the files can be compiled from markdown to html, pdf etc.. using a tool like [pandoc](http://pandoc.org/).  `Rstudio` also enables this automatically. I will do a demonstration in class. The downside of markdown (in comparison to the more complicated markdown languages like html or LaTeX) is that the level of sophistication for formatting of the documents is quite limited (but is being developed). It is usually good enough for an electronic lab notebook, or reports generated to share with your lab members, or as a supplement for a paper to share all of your code and analysis. However, markdown can except both html and LaTeX formatting (and can convert markdown to these formats) and as such can be used to produce beautiful documents.

## R flavoured markdown
R flavoured markdown (`.Rmd`) extends the idea of markdown by enabling the running of `R` code, generating analyses, tables, plots, etc... that can be integrated into your markdown document. While I personally do all of this using a combination of the library `knitr` using the `knit()` function followed by using `pandoc`, `Rstudio` enables all of this easily enough as I will demonstrate in class.

## A few words on project organization

One of the most important things you can do from the get go is organize your project. What you want to avoid is this:

![alt text](http://www.phdcomics.com/comics/archive/phd052810s.gif "PhD Comics 1323, Copyright Jorge Cham")

So the best thing to do is decide from the outset on a useful and clear way to organize your project and your files (and using version control to avoid the problem above). Your book (Bioinformatics Data Skills) goes over this in detail, and other good resources are [here](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424), [here](http://nicercode.github.io/blog/2013-04-05-projects/), or [here](http://www.carlboettiger.info/2012/05/06/research-workflow.html). However, here is my own short and sweet version of it. Each project should be organized into its own folder, with a specific (and common for all projects) set of sub-folders:

```bash
/ProjectName
  ./data
  ./scripts
  ./outputs
  ./misc
  ./manuscript
  README
```

- The `data` folder contains the raw data (read-only!!), hopefully as flat text files (i.e. ``.csv`, .fastq). You should  not alter the data directly here, but instead any editing of the data should be done in the scripts (where you can also explain any changes made to your data). It is ok to have additional sub-folders here to help organize your files.
- The `scripts` folder contains each of your analysis scripts. Instead of writing one very long script, sometimes it is best to break it apart into smaller scripts. One to munge and clean data. For one data quality control (and looking for outliers) and one for the analyses themselves. If you have to create a large number of custom functions, it is good to have a seperate source file. Some computational biologists also prefer to have one additional folder in the hierarchy for these (`./src`).
- The `outputs` folder contains outputs (figures, tables, reports) that can be automatically regenerated from the scripts and the data. THis is not the place for figures that need to be manually (i.e. inkscape, Illustrator) worked on.
- The `misc` folder is for miscellaneous files (like figures you need to manually fix) that can not be automatically regenerated. It is also potentially a good spot to keep your e-lab notebook for the project (some do it as part of the README)
- Finally, the `README` is a file that contains all of your basic meta-data including information on the provenance of the data itself (including experimental design, who collected it, where and when) and other important information (if the data was from another source, links to the DOI, etc..)``

## Writing your code with *style*

For most programming languages there is considerable flexibility on how to you might write your code. However, that can easily mean that everyone uses very different philosophies for naming of variables, functions, etc... Far worse is when an individual changes their approach to naming and coding conventions within a script. This can make even some of the most readable languages (like `python`) very unreadable. So many use *syntax style guides* for programming. This helps to keep your stylistic choices consistent and makes your scripts much more readable. Not only for other people, but for you 6 months down the road when you have completely forgotten what you did. Probably the main two contenders for a syntax style guide are from Hadley Wickham (who develops many libraries in R), which can be found [here](http://adv-r.had.co.nz/Style.html) and the google R style guide [here](https://google.github.io/styleguide/Rguide.xml). I also have my own style guide that I have used in the past for teaching, based largely on the google style guide which is [here](https://msu.edu/~idworkin/ZOL851_style_guide.html). The choice is yours, but use one, and **BE CONSISTENT**.


## Version control using git

First off if you have a copy of Vince Buffalo's book [Bioinformatics Data Skills](http://shop.oreilly.com/product/0636920030157.do), read the chapter "Git For Scientists". (Also read chapter 2 about organizing projects.. Just read the whole book really!!!).

Version control is about how to manage and keep track of editing and revisions on any kind of document. This can be code for software, your analysis pipeline, or manuscript editing. The basic idea is that without keeping track of edits and revisions you make can be problematic.

Imagine a situation where you have some functions in a file you wrote `my_functions.R` you have written that work, but are very slow and do not scale up to the genomic analysis you are planning on making. So you revise the computer code to be more efficient (code refactoring), but do the same thing. So you overwrite your
 original file with the new functions. You now use the functions in this file and see that during your revisions of the functions, you have introduced some unintended behaviour in it. So you want to go back to your original function, but unfortunately you have already overwritten it, so you need to start from scratch.

So next time you are smarter, and when you revise your code you make a new file `my_functions_version2.R`. Everything seems hunky dory, but you realize that some functions are still too slow, and you have learned some new programming techniques to speed things up. As some of these programming techniques are new to you, you are very careful about saving all of your changes in files. However you now end up with a series of files `my_functions_version2.R` to `my_functions_version10.R`. Sadly you realize that one of your new functions has a bug in it, and does not produce the correct results. unfortunately you do not know where you introduced the bug, so it could be anywhere between `my_functions_version3.R` and `my_functions_version10.R`. You end up spending a lot of time tracking this down.

The point of version control is to manage all of these issues by keeping track of the changes you make to a file (the differences), when you made them and also to force you to make a small comment associated with each change to help you find when you made a particular change. This also enables collaborations with other people (with respect to revising code or documents), as you can see what changes they made, and can choose to accept them (or not).
  There are many version control systems out there, but git is a very popular one, and as we will see it is widely used in genomics as we can easily manage projects with online systems like [github](https://github.com/). Indeed github has a great aspect to it that it renders markdown documents (.md) files, enabling us to get some basic formatting (which is what we have been using for this class.)

## learning git and github.

Now I am at most a basic user of both git and github. Everything I have learned is from the tutorials that they provide (and it is generally enough to get started).

First we need to [setup git](https://help.github.com/articles/set-up-git/)

Next we need to [create a repo](https://help.github.com/articles/create-a-repo/). Importantly the name of the repo must match exactly the name of your folder on your local machine as well.

We will now go ahead and make a local repository (*repo*) at the command line.

We can actually follow along with the instructions on github, by clicking on creating a new repo. Give it a name. For the moment **DO NOT** initialize the repo with a README.

It will actually give you the commands you need to copy and paste.


But first at the command line create a folder

```{bash newdir}
mkdir MyNewProject1
cd MyNewProject1
```

Then you can copy and paste what they have.

```{bash copy_git_cmd}
echo "# MyNewProject1" >> README.md
git init
git add README.md
git commit -m "first commit"
git remote add origin https://github.com/DworkinLab/MyNewProject1.git
git push -u origin master
```

the `echo` command is just adding a bit of text to a new file README.md. We will edit this after.
`git init` is initializing a new repo on your local machine.
`git add README.md` is telling the computer to add the file README.md to the repo.
`git commit -m "first commit" is telling the computer to commit changes.
`git remote add origin https://github.com/DworkinLab/MyNewProject1.git` is linking the location of your local repo to github.
`git push -u origin master` is pushing the changes we have made locally to github. I always do all of my changes locally, and push to github to avoid issues. When you are working on collaborations, that will not be possible (see your text)

If you go to your github page, you will now see the changes you have made in the repo.

## making changes.

Let's go ahead and make some local changes. In this case we want to edit the `README.md` file.

Assuming your command line is still in the same directory, go ahead and open the file with your texteditor.

```{bash}
open -a atom README.md
```

Now make some changes to the file.

Something like this:

```{md}
# My first repo
(one number sign acts as a main header, two as a sub-header...)

This is my readme for my first git repo. I am *totally excited* by this. I am also **really** glad Ian likes to write tutorials!.

(put in some code, using md formatting)

That's it
```
Once you have finished that, save the changes and close the file.

Now at your command line type:

```{bash}
touch newfile
```

Which makes a new empty file (just for demonstration purposes, you don't need to do this each time!)


We are now ready to commit our changes:

```{bash}
git add .
git commit -m "edited README, added newfile"
git push
```

Now go ahead and check things out on the github repo!

## **IGNORE BELOW FOR NOW **

## Markdown itself

## bits to incorporate
## Getting started using git and guthub
- What you [don't want](http://www.phdcomics.com/comics/archive.php?comicid=382)
- version control for scripts (and small data)
-

## Reproducible research
- markdown vs. LaTeX
-- [mardown tutorial](http://daringfireball.net/projects/markdown/)
-- [github flavoured markdown](https://help.github.com/articles/github-flavored-markdown)
-- [github flavoured markdown cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)
- knitR & sweave

## Scripting
- Philosophy for organizing scripts.
- Avoid doing [this](http://www.phdcomics.com/comics/archive.php?comicid=1323).
- Source scripts (for functions) and analysis scripts.
- Syntax style guide for [R](https://www.msu.edu/~idworkin/ZOL851_style_guide.html)
- Syntax style guide for python

## Some brief notes on *Sanity Checks* during analysis
- Philosophy: Assume there are mistakes in the data and analysis until you convince yourself otherwise.
- Sanity checks on the computational process (Unit testing)
- Sanity checks on data (labeling, units, extra zeroes...)
- Sanity checks on the analysis
- Some important readings.
