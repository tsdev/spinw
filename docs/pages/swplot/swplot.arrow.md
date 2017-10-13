---
{title: swplot.arrow, link: swplot.arrow, summary: creates a 3D arrow patch, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_arrow, folder: swplot, mathjax: 'true'}

---
  
### Syntax
  
`hPatch = swplot.arrow(rStart, rEnd, R, alpha, lHead)`
 
`hPatch = swplot.arrow(rStart, rEnd, R, alpha, lHead, nPatch)`
  
`hPatch = swplot.arrow(handle, ...)`
 
### Description
  
`hPatch = swplot.arrow(rStart, rEnd, R, alpha, lHead)` draws 3D arrows
between a given start and end position. The arrows will be a triangulated
[patch](https://www.mathworks.com/help/matlab/ref/patch.html) object.
   
`hPatch = swplot.arrow(rStart, rEnd, R, alpha, lHead, nPatch)` creates 
arrows with $$5 n_{patch}$$ number of patch faces per arrow.
 
`hPatch = swplot.arrow(handle, ...)` adds the generated patch object to a
given axis if `handle` is an axis handle or adds the arrows to an
existing [patch](https://www.mathworks.com/help/matlab/ref/patch.html) object, if the given `handle` points to a patch
object.
   
### Examples
 
Draw a 100 random arrows in the $$(\pm 1,\pm 1,\pm 1)$$ cube:
 
```matlab
swplot.figure
N = 100
swplot.arrow(2*rand(3,N)-1,2*rand(3,N)-1,0.01,30,0.05)
swplot.zoom(40)
```
 
{% include image.html file="generated/swplot__1.png" alt="swplot.zoom(40)" %}
 
### Input Arguments
  
`handle`
: Handle of an axis or triangulated patch object. In case of patch
  object, the constructed faces will be added to the existing object.
  
`rStart`
: Coordinates of the arrow starting point, one vector per arrow in a
  matrix with dimensions of $$[3\times n_{obj}]$$.
  
`rEnd`
: Coordinates of the arrow end point, one vector per arrow in a
  matrix with dimensions of $$[3\times n_{obj}]$$.
  
`R`
: Radius of the arrow body, scalar.
  
`alpha`
: Angle of the head in degree.
  
`lHead`
: Length of the head.
  
`nPatch`
: Number of points on the circle of the body, default value is stored in
  `swpref.getpref('npatch')`. The final patch object will have
  $$5n_{patch}$$ number of faces and $$3n_{patch}$$ number of vertices.
  
### See Also
  
[swplot.cylinder](swplot_cylinder)
 

{% include links.html %}
