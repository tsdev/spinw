---
{title: swplot.cylinder( ), summary: draws a closed/open 3D cylinder, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_cylinder.html, folder: +swplot, mathjax: 'true'}

---
draws a closed/open 3D cylinder
 
hPatch = SWPLOT.CYLINDER(rStart, rEnd, R, {nPatch}, {close})
 
The function can draw multiple cylinders with a single patch command for
speedup.
 
hPatch = SWPLOT.CYLINDER(handle,...)
 
Handle can be the handle of an axes object or a patch object. It either
selects an axis to plot or a patch object (triangulated) to add vertices
and faces.
 
Input:
 
handle    Handle of an axis or patch object. In case of patch object, the
          constructed faces will be added to the existing object instead
          of creating a new one.
rStart    Coordinate of the starting point with dimensions [3 nCylinder].
rEnd      Coordinate of the end point with dimensions [3 nCylinder].
R         Radius of the arrow body.
nPatch    Number of points on the curve, default value is stored in
          swpref.getpref('npatch').
close     If true the cylinder is closed. Default is true.
 
See also SWPLOT.ARROW.
 
