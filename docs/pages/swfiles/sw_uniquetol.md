---
{title: sw_uniquetol, link: sw_uniquetol, summary: returns the unique column vectors
    within tolerance, keywords: sample, sidebar: sw_sidebar, permalink: sw_uniquetol,
  folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`[Mu, firstIdx] = sw_uniquetol(M,tol)`
  
### Description
  
`[Mu, firstIdx] = sw_uniquetol(m,tol)` returns unique column vectors
within the given `tol` tolerance. Two column vectors are considered
unequal, if the distance between them is larger than the tolerance
($$\delta$$):
 
$$\sqrt{\sum_i (V_i-U_i)^2} < \delta$$
  
### Input Arguments
  
`M`
: Matrix that contains column vectors.
  
`tol`
: Distance tolerance, default value is $$10^{-5}$$.
  
### Output Arguments
  
`Mu`
: Matrix that contains the unique column vectors.
 
`firstIdx`
: Indices pointing to the first occurence of the unique element.
 
This function is similar to the Matlab built-in
`unique(M,'rows','first')`, but with controllable tolerance.
 
### See Also
 
[unique](https://ch.mathworks.com/help/matlab/ref/unique.html)
 

{% include links.html %}
