---
{title: sw_nvect, link: sw_nvect, summary: determines the best normal vector for the
    set of vectors, keywords: sample, sidebar: sw_sidebar, permalink: sw_nvect, folder: swfiles,
  mathjax: 'true'}

---
  
### Syntax
  
`[n, collinear] = sw_nvect(V)`
  
`[n, collinear] = sw_nvect(V,epsilon)`
 
### Description
  
`[n, collinear] = sw_nvect(V)` determines whether the given set of
vectors are collinear or coplanar. If they are coplanar, it returns the
best fitting normal vector, while if they are collinear returns the
average of the given vector.
 
The function can also deal with complex vectors, separating the real and
complex parts as separate vectors.
 
 
`[n, collinear] = sw_nvect(V,epsilon)` also gives the upper limit of the
collinearity controlled by `epsilon`.
   
### Input Arguments
 
`V`
: Matrix of column vectors with dimensions of $$[3\times N]$$. Where each
  column defines a vector.
 
`epsilon`
: Defines the limits of collinearity with the following values:
  * `1`   the function always return the `n` closest to the collinear
          direction,
  * `2`   the function always return the `n` vector closest to the normal 
          of the coplanar plane.
  * `e`   upper limit of collinearity, default value is 0.1, smaller
          positive values mean stricter limits on collinearity.
   
### Output Arguments
 
`n`
: Row vector parallel to the collinear vector direction or
  perpendicular to the best fitting plane of the coplanar vectors.
   
`collinear`
: `true` if the given set of vectors are collinear.
   

{% include links.html %}
