---
{title: sw_mirror, link: sw_mirror, summary: mirrors a 3D vector, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_mirror, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`[~, M] = sw_mirror(n)`
  
`[Vp, M] = sw_mirror(n,V)`
 
### Description
  
[~, M] = sw_mirror(n) generates the transformation matrix corresponding
to a mirror plane perpendicular to `n`.
 
`[Vp, M] = sw_mirror(n,V)` mirrors the vectors in `V`.  
  
To mirror any column vector use the following:
 
```matlab
Vp = M * V
```
 
To apply mirror plane operation on tensors ($$3\times 3$$ matrices) use the
following command:
 
```matlab
Ap = M * A * M'
```
 
### Input Arguments
  
`n`
: 3D Row vector, normal to the mirror plane.
  
`V`
: Matrix of 3D vectors, dimensions are $$[3\times N]$$.
  
### Output Arguments
  
`Vp`
: Mirrored vectors in a matrix with dimensions of $$[3\times N]$$.
 
`mirM`
: Matrix of the mirror transformation, dimensions are $$[3\times 3]$$.
  
### See Also
  
[sw_rot](sw_rot)
 

{% include links.html %}
