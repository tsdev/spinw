---
{title: swplot.circle, link: swplot.circle, summary: creates a circle surface in 3
    dimensions, keywords: sample, sidebar: sw_sidebar, permalink: swplot_circle.html,
  folder: swplot, mathjax: 'true'}

---
 
hPatch = SWPLOIT.CIRCLE(r0, n, R, {N})
 
hPatch = SWPLOIT.CIRCLE(handle,...)
 
Handle can be the handle of an axes object or a patch object. It either
selects an axis to plot or a patch object (triangulated) to add vertices
and faces.
 
Input:
 
handle    Handle of an axis or patch object. In case of patch object, the
          constructed faces will be added to the existing object instead
          of creating a new one.
r0        Center of the circle, vector with three elements.
n         Vector normal to the circle surface, vector with three elements.
R         Radius of the circle.
N         Number of points on the curve, default value is stored in
          swpref.getpref('npatch').
 
Example:
 
swplot.circle(zeros(3),eye(3),1,100)
 
See also SWPLOT.CYLINDER.
 

