# R exercise 6

0. If you are not working with a fresh R instance, then use the command:
```R
rm(list=ls()
```
1. Create a variable `a` for a numeric vector with the following values:
  - ` 3 9 12 17 8 24 3 9 34`
2. Compute the coefficient of variation (standard deviation divided by the mean) for `a`.
3. Now make a function called `CoefVar` that enables you to calculate the coefficient of variation for any arbitrary numeric vector. Try it on `a` and confirm you get the same value for the coefficient of variation.
4.  Now use your function `CoefVar` on the following vector:
```R
b <- c(-4, -2, 0, 1, 3, -5, -1, -8)
```
5. What is your result? The result likely does not make sense (check the wikipedia [entry](https://en.wikipedia.org/wiki/Coefficient_of_variation)). 
6. Fix your function `CoefVar` to fix the bug in the calculation.
7.  Save your script you just wrote (with a name like `MyFirstFunction.R`)
8. Close `R` (do not save your workspace).
9. Re-open `R` and source in your file `MyFirstFunction.R`. Be careful about making sure it is in the same working directory as you are in using `getwd()`, or setting it with `setwd()`.
10. How can you check to see if the function `CoefVar` is in your environment (your workspace)? 
11. Run the function on a numeric vector to make sure it works.