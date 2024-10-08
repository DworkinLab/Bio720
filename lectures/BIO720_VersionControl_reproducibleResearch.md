---
title: "Reproducibel research, part I"
author: ""
date: "03 Sep 2024"
output:
  slidy_presentation: 
    keep_md: yes
    highlight: tango
    fig_retina: 1
  ioslides_presentation: default
editor_options: 
  chunk_output_type: console
---



# Reproducible research is key to better science!

## link to github page for the class on this
[go here](https://github.com/DworkinLab/Bio720/blob/master/IntroductionMarkdownAndVersionControl/Bio720_IntroductionMarkdown.md)

## Learning Objectives for this week


## Readings and DataCamp
- DataCamp, I will assign the courses on RMarkdown and git. ON DataCamp 

- Read [Good enough practices in scientific computing](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005510) by Wilson et al.  

- Read chapters 1 and 2 in [Bioinformatics Data Skills](https://mcmaster.primo.exlibrisgroup.com/permalink/01OCUL_MU/deno1h/alma991030055129707371)

- An alternative reading on version control and using git and github is [chapter 2 in the CSB] book(https://doi-org.libaccess.lib.mcmaster.ca/10.1515/9780691183961-006). Pick which one you prefer.


## What's wrong with this?

![alt text](http://www.phdcomics.com/comics/archive/phd052810s.gif "PhD Comics 1323, Copyright Jorge Cham")

## How do we make science more reproducible?

- Is publishing your paper enough?

## How do we make science a more open and reproducible endeavor?
- Data?  
- Meta-data?   
- Analysis scripts/code?

What easy steps can make science more reproducible & helpful for future scientists?

## How do we make science a more open and reproducible endeavor?
- You want to make this process as easy as possible for yourself first and foremost.
- In the process of organizing your data and analysis pipelines to share, you will make yourself a better scientist.
- Future you (6 months after you finish data collecting), will thank you too!


## Tri-council (NSERC, CIHR and SSHRC) guidelines are in place
"Grant recipients are required to deposit into a recognized digital repository all digital research data, metadata and code that directly support the research… 
The repository will ensure safe storage, preservation, and curation of the data. The agencies encourage researchers to provide access to the data ... 
Whenever possible, these data, metadata and code should be linked to the publication with a persistent"
[See links here](https://www.ic.gc.ca/eic/site/063.nsf/eng/h_97610.html
)


## So how do you make this easy on yourself?

- At the onset of your research, organize your workflow with the goal of making it as easy as possible to share.

- With a few simple guidelines, you will save yourself an enormous amount of time while making your data and your analysis scripts available.

## Rule 1:Organize your project in a simple and clear manner
```bash
/ProjectName
  ./data
  ./scripts
  ./outputs
  ./misc
  ./manuscript
  README
```
## Rule 1:Organize your project in a simple and clear manner

All projects should have very similar organizational formats:

- Informative project name.  
- `data` folder contains raw data as a flat file (tabular data, comma seperated), database or if necessary a proprietary spreadsheet.  

- `scripts` folder contains scripts/code you use to analyse your data.

- `outputs` for automated outputs (tables and figures from the scripts + data)  

- `misc` for miscellaneous. Additional meta-data for files, experimental design information, readme of overall project. Figures that can not easily be generated.

## Rule 2: strive to store raw data as plain text files
- File formats for many commercial programs (SAS, JMP, Excel,) can not be easily read by other programs.

- Worse still. As new versions of software emerge, no guarantee of being able to open older file formats (Excel).
-- It is fine to use Excel. But when sharing data (and scripts) use formats like .txt (such as tab delimited or .tsv) or .csv (comma seperated) which can be easily read.

## Rule 2.5: Consistent file naming of raw data
- Many projects generate hundreds or thousands of files.
For instance digital microscopy, automated measuring tools, sequencing, etc.
- Spend a few minutes thinking about how to name your files. 
This will save you hours or days down the road!

## Rule 2.5: Consistent file naming of raw data
Some "do's" and "do not"

**DO**: Use the same naming conventions for every file. i.e.

ID_Dmel_M_HN_01.tif   
ID_Dmel_M_HN_02.tif  
ID_Dmel_M_HN_03.tif  
ID_Dmel_F_HN_01.tif  
ID_Dmel_F_HN_02.tif  
ID_Dmel_F_HN_03.tif 

## Rule 2.5: Consistent file naming of raw data
 **DO NOT**:Use spaces in file names. Use consistent delimiters, and limit it to one type of delimiter if you can. i.e. underscores!

ID_Dmel_M_HN_01.tif  

**Not** ID Dmel_M_HN 01.tif

**Not** ID_Dmel-M_HN.01.tif 

## Rule 2.5: Consistent file naming of raw data

**DO**: Use consistent delimiters! Underscore `_` and dash `-` are good choices. 

## Rule 2.5: Consistent file naming of raw data

**DO**: Try to have as much of the experimental information in the file name (like below)

ID_Dmel_M_HN_01.tif 

ID = initials  
Dmel = species name  
M = sex (M/F)  
HN = High Nutrition  
01 = specimen number  

## Rule 2.5: Consistent file naming of raw data

Why is this important?

It allows for all of the experimental variables to be extracted in an automated way (even in Excel).

But having inconsistencies in naming can make a 1 minute activity or one line of code into something not fun.

At least 40% of my coding time is spent parsing and cleaning badly organized data, in particular due to naming conventions.


## Rule 3: version control analysis (and if possible, data)
- Version control is about how to manage and keep track of editing and revisions on any kind of document. This can be code for software, your analysis pipeline, or manuscript editing. 

- The basic idea is that without keeping track of edits and revisions you make can be problematic.

## version control analysis (and if possible, data)
- version control keeps track of changes you make to files (differences), when you made them. It also forces you to make a small comment associated with each change to help you find when you made a particular change. 

- This also enables collaborations with other people (with respect to revising code or documents), as you can see what changes they made, and can choose to accept them (or not). 

- There are many version control systems out there. 

## Rule 4 (or rule 1!) Write a readme file right at the beginning of your analysis
- Don’t assume that any collaborators or future users (or future you) will know what variables are in your data set.

- Write a short readme file containing basic meta-data for the project including: explanations of the variables, how and where they were collected (and why), file naming conventions etc.

- If you do this while organizing everything it takes just a few minutes, but it will save future you and collaborators so much time.

## Sharing data (and code)
In addition to backing up data for yourself, getting it out there as soon as possible (in conjunction with manuscript submission or acceptance) is crucial. 


## Where to share data
- There are often both institutional resources for long term data archiving, and broader initiatives.

- Some of the Data sharing portals are fairly general purpose (DRYAD, zenodo, figshare, dataverse, FRDR in Canada).

- Some are highly specialized for certain data types (Morphobank, tree of life, NCBI GEO,NCBI SRA, morphosource, DigiMorph

- [A useful short blog post about this](https://evodify.com/free-research-repository/)













