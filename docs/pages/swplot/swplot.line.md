---
{title: swplot.line, link: swplot.line, summary: draws a 3D line using patch, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_line.html, folder: swplot, mathjax: 'true'}

---
 
hLine = SWPLOT.LINE(rStart, rEnd, {lineStyle}, {lineWidth},{multiPatch})
 
Plots line disconnected segments between multiple rStart --> rEnd pairs
of coordinates.
 
hLine = SWPLOT.LINE(r, [], {lineStyle}, {lineWidth},{multiPatch})
 
Plots multiple continuous curves defined in the r matrix.
 
hLine = SWPLOT.LINE(handle,...)
 
Handle can be the handle of an axes object or a line object. It either
selects an axis to plot or a patch object (triangulated) to add vertices
and faces.
 
Input:
 
handle    Handle of an axis or patch object. In case of patch object, the
          constructed faces will be added to the existing object instead
          of creating a new one.
rStart    Coordinate(s) of the starting point, either a 3 element vector or
          a matrix with dimensions [3 nLineSegment] to plot multiple line
          segments.
rEnd      Coordinate(s) of the end point, either a 3 element vector or
          a matrix with dimensions [3 nLineSegment] to plot multiple line
          segments.
r         Matrix with dimensions [3 nCurve nPointPerCurve]. The function
          will plot a nCurve number of disconnected curves. The i-th
          curve will follow the x=r(1,i,:), y=r(2,i,:), z=r(3,i,:)
          (parameteric) curve.
 
lineStyle Line style, default is '-' for continuous line.
lineWidth Line with in pt, default is 0.5.
mPatch    If true, a separate patch object will be created per line
          segment. Default is false, a single patch object will store all
          line segments.
 
See also LINE.
 

