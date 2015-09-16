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

 We also learned to navigate around to other directories (i.e. `cd /` to change to the root directory, or `cd ~` to go back to the home directory). We discussed making new files (`touch newFile`), removing files (`rm newFile`), making new directories (`mkdir myNewDirectoryName`) and removing directories (`rmdir`). Finally we used `more` to look at the content of a file one page or line at a time.
