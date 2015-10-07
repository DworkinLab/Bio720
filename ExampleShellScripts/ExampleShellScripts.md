#Example Shell Scripts (Bio720)

In class today (October 6th 2015) we talked a bit how the three different instructors all use *shell scripts* to run our computational jobs. To that end, we wanted to provide you with a few examples shell scripts for some of the common steps you will use at the command prompt. For a bit more of an introduction to shell scripts,  consult either the *Practical computing for biology" book or look at the online tutorial at software carpentry on [shell scripts](http://swcarpentry.github.io/shell-novice/05-script.html).

## What is a shell script?

A *shell script* is simply a collection of unix commands organized into a file. When done this way, the whole script can be run, instead of having to do everything interactively (which is useful at the beginning, but a pain and inefficient once you have a working pipeline). The only thing seperating a shell script from a regular set of commands is that you have to let your shell know it is a shell script. You do this in two ways. First (although this is mostly for human readability) we generally give shell scripts the `.sh` extension (for *shell*).

Second, and most importantly, we let the computer know which *shell* we want to use (in class we have exclusively been using *bash*), so we tell the computer the location of the bash shell program. Usually that means the first line of a shell script looks like:

```bash
#!/bin/bash
```

Which tells the computer to look for the bash program in the `/bin` folder (for binaries). After that each line is just the standard commands.


## An example with trimmomatic

Today in class BG went over how you could run [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), which is one of a number of tools that can be used to remove both adapter sequences and low quality sequences from your `.fastq` files. As we (the instructors) all use slightly different computational resources to run our scripts, we each use slightly different shell scripts to run this.

For instance BE uses a script like [this](./BE_trimmomatic.sh). This is a pretty standard approach to running many files through the same program (or pipeline of programs in some cases). I will not go through each line here (as we have discussed some of this in class, and the rest through our suggested readings and tutorials). However a few salient points.

```bash
files=”Heart
Liver
Lung
Muscle
Kidney
Ovary”
```

Here BE is generating a single variable (`files`) which is storing the names of all of sequence files (which he has previously renamed Heart, Liver,....Ovary). This allows him to loop through (using the `for` loop in the next part of the script) each set of paired files (these are paired end reads as you can tell by the `PE` used in the program) and run trimmomatic. For more information on using the `for` loop in bash, see [this tutorial](http://swcarpentry.github.io/shell-novice/04-loop.html).

In the `for` loop, you can see that instead of providing the name of each sequence file (Heart, Liver, etc) he instead just writes `${file}`. The `for` loop is going through each element (in this case the names of each of the files as named by tissue type) and substitures `${file}` with the actual file name. This is an example of using a *variable* (which is denoted with the `$`). This is a standard approach in most programming languages.

One thing worth noting is that when running the shell script above, you could try something like:
```bash
./BE_trimmomatic.sh
```

However, you may (depending on the system) get an error. This is because the `.sh` means that the computer is running this as its own program. However, with UNIX you may need to make sure you have your permissions set to run this program. I would definitely say for your own scripts, this is probably a useful template.

## Working with a scheduling system

As an alternative, here is an [example of a script](sd_BSA2014_trimmomatic.sh) I would have submitted to run trimmomatic. In my case, the shell script looks a bit different (in particular the first few lines) because I have submitted it to a large high performance computing cluster with ~6000 cores, across several hundred nodes (and many slightly different architectures). Hundreds of researchers use this system, and as such they have to develop some system to schedule and prioritize submitted jobs. Thus the first few lines starting with the `#` (which would normally represent a commented out line) provide information to the scheduling program to state how many resources I need, and for how long. The particulars of that do not matter so much. The one nice aspect of using a scheduler like that is it is easier to *spawn* many serial jobs concurrently (i.e. instead of running each job sequentially, I can run many jobs at the same time, depending on the constraints of the computational resources). As an alternative to using a `for` loop, I have used an array to spawn these jobs. With the PBS scheduler, I set one flag

```bash
#PBS -t 0-7
```
Which let's it know I want to run 8 jobs (called 0 1 2 ... 7)

In the script itself I set up an *array* variable here (using a wild-card to grab all files ending in *.fastq.gz):
```bash
files=(${dir}*_R1.fastq.gz)
```

and then can pull out individual elements as I establish the formal array
```bash
file=${files[${PBS_ARRAYID}]}
```
with the `PBS_ARRAYID` just being the numbers from 0-7.

I then go through a few little tricks to simplify the names of each file and folder, simply to make the trimmomatic call a bit cleaner to read.
