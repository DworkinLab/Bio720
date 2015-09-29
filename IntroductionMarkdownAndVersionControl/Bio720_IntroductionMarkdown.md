# Getting used to markdown

## What is markdown?

Markdown is a very simple and readable file format that enables the generation of *dynamic documents* that allow for the easy generation of reports, presentations and even drafts of papers. Importantly it makes embedding blocks of code really easy, so it can help you if you want to keep an electronic copy of your lab notebook. As we will see, with version control (like `git`) this can enable a powerful way of organizing your projects.  Since markdown is so simple, it really only takes about 5 minutes to get started, and you will be comfortable with it after a couple of days. There are many introductions to markdown on the web. I usually recommend [this one](https://daringfireball.net/projects/markdown/basics), although if you decide (after I introduce it) to use github, I suggest looking at [this as well](https://help.github.com/articles/markdown-basics/).

The goal of `markdown` is to keep everything human readable even as it allows for some simple formatting and generation of figures that can be imbededded. While the markdown (using the extension `.md` files are just regular text files, they can be converted to PDF, html and many other formats easily. Indeed a number of text editors (google "markdown editor") can render the markdown in realtime. Otherwise the files can be compiled from markdown to html, pdf etc.. using a tool like [pandoc](http://pandoc.org/).  `Rstudio` also enables this automatically. I will do a demonstration in class.

## R flavoured markdown
R flavoured markdown (`.Rmd`) extends the idea of markdown by enabling the running of `R` code, generating analyses, tables, plots, etc... that can be integrated into your markdown document. While I personally do all of this using a combination of the library `knitr` using the `knit()` function followed by using `pandoc`, `Rstudio` enables all of this easily enough as I will demonstrate in class.

## A few words on project organization

One of the most important things you can do from the get go is organize your project. What you want to avoid is this:

![alt text](http://www.phdcomics.com/comics.php?f=1323 "PhD Comics 1323, Copyright Jorge Cham")

So the best thing to do is decide from the outset on a useful and clear way to organize your project and your files (and using version control to avoid the problem above). Your book (Bioinformatics Data Skills) goes over this in detail, and another good resource is [here](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424). However, here is my own short and sweet version of it. Each project should be organized into its own folder, with a specific (and common for all projects) set of sub-folders:

```bash
ProjectName
  ./data
  ./scripts
  ./outputs
  ./misc
  ./manuscript
  README
```

- The `data` folder contains the raw data (read-only!!), hopefully as flat text files (i.e. ``.csv`, .fastq). You should endevour to not alter the data directly here, but instead any editing of the data should be done in the scripts. It is ok to have additional sub-folders here to help organize your files.
- The `scripts` folder contains each of your analysis scripts. Instead of writing one very long script, sometimes it is best to break it apart into smaller scripts. One to munge and clean data. For one data quality control (and looking for outliers) and one for the analyses themselves. If you have to create a large number of custom functions, it is good to have a seperate source file. Some computational biologists also prefer to have one additional folder in the hierarchy for these (`./src`).
- The `outputs` folder contains outputs (figures, tables, reports) that can be automatically regenerated from the scripts and the data. THis is not the place for figures that need to be manually (i.e. inkscape, Illustrator) worked on.
- The `misc` folder is for miscellaneous files (like figures you need to manually fix) that can not be automatically regenerated. It is also potentially a good spot to keep your e-lab notebook for the project (some do it as part of the README)
- Finally, the `README` is a file that contains all of your basic meta-data including information on the provenance of the data itself (including experimental design, who collected it, where and when) and other important information (if the data was from another source, links to the DOI, etc..)

## Markdown itself
