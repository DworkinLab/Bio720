# Logging onto a remote machine

In this tutorial we are going to log on to a remote machine at the terminal. Most computational work you will do requires more processing power and memory than on your local work station (i.e. your laptop). For this tutorial we will use the Golding lab cluster (Thanks Brian!).

## Open your terminal.

If you are using either a computer running a linux operating system or Mac OS X, you are already using an operating system with Unix. 

### For a Mac

If you do not know where your terminal already is (Applications > Utilities), you can use spotlight to find it. Just starting typing `terminal` in spotlight and it will show up. You may want to put the terminal application in your application shortcuts in the dock (you will use it a lot in this class).

### For Linux

If you are already using Linux, chances are you can already find your terminal!

# For Windows

Windows does not natively run any flavour of UNIX. You have a copy of easy options. The classic tool is [`putty`](http://www.chiark.greenend.org.uk/~sgtatham/putty/) which allows you to log onto remote machines and use the terminal. There are also UNIX terminal emulators such as [CygWin](https://www.cygwin.com/). FInally, you can also make your local machine dual boot (so that you can boot up either Windows or Linux), but that is outside the pervue of the course (if you do try it, back up your data first).

## Connecting to a remote machine using `ssh`

Now we are going to connect to the Golding Lab cluster via the terminal. We are going to use a simple program within UNIX called `ssh` which stands for SecureShell. 

Hopefully you have already received your login credentials for the Golding Lab cluster.  If you have then you want to type out the following line at the command prompt:

```
ssh yourUserName@info.mcmaster.ca
```

Where you will substitute your own username for YourUserName.

You should then see a line that looks like the following:

```
yourUserName@info.mcmaster.ca's password:
```

where you type in your password (and then press return). You will get three tries after which you will be locked out for 20 minutes.

Once you successfully login, you should see something like the following:
```
Last login: Mon Aug 17 09:55:50 2015 from 13-15-64.client.wireless.msu.edu
[yourUserName@infoserv ~]$ 
```
Where the `$` represents the command prompt, and you can use it like you are on any other system.

## Navigating the folder structure
