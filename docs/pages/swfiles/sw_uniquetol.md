---
{title: sw_uniquetol( ), summary: returns the unique column vectors within tolerance,
  keywords: sample, sidebar: sw_sidebar, permalink: sw_uniquetol.html, folder: swfiles,
  mathjax: 'true'}

---
returns the unique column vectors within tolerance
 
[unique, firstIdx] = SW_UNIQUETOL(M,tol)
 
Two column vectors are considered unequal, if the distance between them
is larger than the tolerance.
 
Input:
 
M         Matrix that contains column vectors.
tol       Distance tolerance, default is 1e-5.
 
Output:
 
unique    Matrix that contains the unique column vectors.
firstIdx  Indices pointing to the first occurence of the unique element.
 
This function is similar to the Matlab built-in unique(M,'rows','first'),
but with arbitrary tolerance.
 
