# 22-Lab2-perflab

Welcome to the much belated repository for `22-Lab2-perflab`!

If you haven't started working on the assignment yet, [you can find the assignment description on the website](https://github.itu.dk/pages/OSC/22-Website/labs/perflab.pdf)

In order for us TA's to provide feedback based on your implementations, we would like if you could hand-in your solution in a similar manner as with [22-Lab1-datalab](https://github.itu.dk/OSC/22-Lab1-datalab).

Just as a refresher, we would like you to:

- Fork this repository (using github.itu.dk's web interface)
- Clone your fork to your local machine (using `git clone`)
- Change to a branch named after your initials (using `git checkout -b initials`)

You have most likely downloaded the assignment elsewhere and already made modifications. In that case, make sure that you are cloning the repository in a different location so that you don't lose your changes.

When your fork has been set up, simply copy your kernels.c file into your fork's folder.
E.g. if you have your old perflab folder in `perflab` and you've cloned to the same directory (and thereby created a new folder `22-Lab2-perflab`) then you can copy the file using the following command (on 
linux):

```
# cp <src file> <dst file>
cp perflab/kernels.c 22-Lab2-perflab
```

By default, the destination path is overwritten when you copy.

Then, you can `cd` into the new folder (`22-Lab2-perflab`) and add the changes to a commit by using 

```
git add kernels.c
```

Finally, you can commit your changes by writing 

```
git commit -m "<some commit message>"
```

and push them to your branch using 

```
git push -u origin initials
``` 

After this, you should be tracking your branch, meaning that you can simply write `git push` to upload your changes.

When you are satisfied with your changes/results you should create a pull request using github.itu.dk's web interface and name it like `initials / attempt #` or something along those lines.

The assignment is due on Sunday night 13 November and we will try to have feedback for you within 1-2 weeks.
