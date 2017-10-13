---
{title: swplot.cylinder, link: swplot.cylinder, summary: creates a closed/open 3D
    cylinder patch, keywords: sample, sidebar: sw_sidebar, permalink: swplot_cylinder,
  folder: swplot, mathjax: 'true'}

---
  
### Syntax
  
`hPatch = swplot.cylinder(rStart, rEnd, R)`
  
`hPatch = swplot.cylinder(rStart, rEnd, R, nPatch, close)`
 
`hPatch = swplot.cylinder(handle, ...)`
 
### Description
  
`hPatch = swplot.cylinder(rStart, rEnd, R)` generates multiple cylinders
with a single triangular patch command. The cylinders are defined by
start and end positions and their radii.
   
`hPatch = swplot.cylinder(rStart, rEnd, R, nPatch, close)` creates 
cylinders with $$4 n_{patch}$$ number of patch faces per arrow.
   
Handle can be the handle of an axes object or a patch object. It either
selects an axis to plot or a patch object (triangulated) to add vertices
and faces.
   
### Examples
 
Draw 100 random cylinders within the $$(\pm 1,\pm 1,\pm 1)$$ cube:
 
```matlab
swplot.figure
N = 100;
swplot.cylinder(2*rand(3,N)-1,2*rand(3,N)-1,0.1,100,true)
swplot.zoom(30)
```
 
{% include image.html file="generated/swplot_cyl_1.png" alt="swplot.zoom(30)" %}
  
### Input Arguments
  
`handle`
: Handle of an axis or patch object. In case of [patch](https://www.mathworks.com/help/matlab/ref/patch.html) object,
  the constructed faces will be added to the existing object instead of
  creating a new one.
  
`rStart`
: Coordinate of the starting point or multiple starting points in a
  matrix with dimensions $$[3\times n_{obj}]$$.
  
`rEnd`
: Coordinate of the end point or multiple end points in a
  matrix with dimensions $$[3\times n_{obj}]$$.
  
`R`
: Radius of the arrow body, scalar.
  
`nPatch`
: Number of points on the circle of the body, default value is stored in
  `swpref.getpref('npatch')`. The final patch object will have
  $$4n_{patch}$$ number of faces and $$2n_{patch}$$ number of vertices.
  
`close`
: If `true` the cylinder is closed at both ends. Default is `true`.
  
### See Also
  
[swplot.arrow](swplot_arrow)
 

{% include links.html %}
