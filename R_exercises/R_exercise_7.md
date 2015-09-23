# R exercise 7

0. Clear your workspace of any lingering objects by using:
```R
rm(list=ls())
```
1. Create a vector `a` of the sequence of numbers of from -100 to 100 (inclusive) by 5s (i.e. -100, -95, ...., 95, 100).
2. What is the length of `a`?
3. Create a vector `b` of the same length as `a` that goes from 1 to `length(a)`
4. Generate the following matrix:
```R
C <- matrix(1:99, nrow=10, ncol=10)
```
5. What happens with `C`? Does `R` still produce an object `C`?
6. How many rows and columns does `C` have?
7. What is the final element (final row by final column) in `C`?  
8. Remove the object `C`
9. Create a new object 'C`:
```R
C <- matrix(1:100, nrow=10, ncol=10)
```
10. How many rows and columns does the new version of `C` have? What is the final element of this new version of `C`?
11. If you have not done so already, think about the difference between these two versions of `C`, and how `R` handles putting a shorter object (i.e. `1:99` into an object with more elements.
12. Extract the complete third row of `C`.
13. Extract the complete fourth column of `C`.
14. Extract the element in the 4th row and 5th column.
15. Replace the element in the 4th row and 5th column with the number `1000`.
16. Extract the elements of `a` where elements are less than 10 **OR** greater than or equal to 86. How many elements do you have in this new vector?
