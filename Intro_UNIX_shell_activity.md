# Bio720: Practical for the introduction to the `UNIX` shell

You will need access to a Unix shell (Max OS X or Linux) or use [MobaxTerm](http://mobaxterm.mobatek.net) to login into the info server.

If you have not already done so, read chapter 3 (remedial UNIX skills) in Bioinformatic Data Skills. Alternatively, you could also use a tutorial like this one on using the shell from [software carpentry](http://swcarpentry.github.io/shell-novice/).

## Review

In class this week Brian introduced some of the basic commands on the UNIX shell to help you login into a remote machine using `ssh`

So first you login into the remote machine:

```{shell}
ssh username@server.ca
```

where the username and server information are on your sheets.

Then we learned some commands to navigate around the file system like:

```{shell}
ls
```

to list the files in the current working directory or:

```{shell}
ls /
```

to look at the contents of the folder.

 We also learned to navigate around to other directories (i.e. `cd /` to change to the root directory, or `cd ~` to go back to the home directory) and determine our working directory (`pwd`). We discussed making new files (`touch newFile`), removing files (`rm newFile`), making new directories (`mkdir myNewDirectoryName`) and removing directories (`rmdir`). We used the `mv` command to both move and rename files, and `cp` to copy files (and `scp` to move files across computers). To get more information about the commands we used (including for finding flags for optional arguments) we used the `man` pages (like `man ls`). we used `more` to look at the content of a file one page or line at a time. To examine the jobs running on the machine we used both `top` and `ps`.
 
## Your activity

Given the tutorial, and your readings, we now want you to do several things.

1. log onto the remote machine (info) that you have access to. What directory do you start in when you login?
2. Without changing your working directory, look at the contents in the root (`/`) directory.
3. the info server has most of the bioinformatic software located in `/usr/local`. Please navigate (change directory) to this location.
4. using what we learned in class (and maybe using `man` to find the right flags) take a look at all of the files. You will notice that there are more items in this directory than can be displayed on the page at one time. Instead look up the arguments for the `ls` command that will allow you to simultaneously do the following:
     display each item in the directory in its own row.
     provide more information about each item (in particular permissions like `drwxr-xr-x`).
5. What is the permission for the directory to the program `bowtie2` for you (located in `usr/local/`).
6. return to your home directory.
7. make a new folder called NewFolderTest.
8. Navigate into NewFolderTest.
9. Create a new empty file called EmptyFile.txt (hint remember the `touch` command).
10. Copy EmptyFile.txt to a new file EmptyFile2.txt.
11. Rename EmptyFile.txt to EmptyFile3.txt
12. navigate back up to your home directory.
13. try to delete NewFolderTest. What happens? How can you delete NewFolderTest and all of its contents?
14. What is the top running process currently on the machine you are working on?

