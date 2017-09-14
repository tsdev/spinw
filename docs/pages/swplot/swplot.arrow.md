---
{title: swplot.arrow, link: swplot.arrow, summary: draws a 3D arrow using patch, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_arrow.html, folder: swplot, mathjax: 'true'}

---
 
hPatch = SWPLOT.ARROW(rStart, rEnd, R, alpha, lHead, {nPatch})
 
hPatch = SWPLOT.ARROW(handle,...)
 
Handle can be the handle of an axes object or a patch object. It either
selects an axis to plot or a patch object (triangulated) to add vertices
and faces.
 
Input:
 
handle    Handle of an axis or patch object. In case of patch object, the
          constructed faces will be added to the existing object instead
          of creating a new one.
rStart    Coordinate of the starting point.
rEnd      Coordinate of the end point.
R         Radius of the arrow body.
alpha     Angle of the head in degrees.
lHead     Length of the head.
nPatch    Number of points on the curve, default value is stored in
          swpref.getpref('npatch').
 
See also SWPLOT.CYLINDER.
 

