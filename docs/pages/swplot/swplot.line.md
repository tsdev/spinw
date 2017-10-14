---
{title: swplot.line, link: swplot.line, summary: creates 3D line patch, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_line, folder: swplot, mathjax: 'true'}

---
  
### Syntax
  
`hLine = swplot.line(rStart, rEnd)`
 
`hLine = swplot.line(r,[])`
  
`hLine = swplot.line(rStart, rEnd, lineStyle, lineWidth, multiPatch)`
 
### Description
 
`hLine = swplot.line(rStart, rEnd)` creates disconnected line segments
between multiple `rStart(:,i)` `rEnd(:,i)` pairs of 3D coordinates. The
lines are shown as patch faces.
   
`hLine = swplot.line(r,[])` creates connected line segments  between
the consicutive points `r(:,i)`.
   
`hPatch = swplot.line(handle, ...)` adds the generated patch object to a
given axis if `handle` is an axis handle or adds the lines to an
existing [patch](https://www.mathworks.com/help/matlab/ref/patch.html) object, if the given `handle` points to a patch
object.  
  
### Input Arguments
  
`handle`
: Handle of an axis or triangulated patch object. In case of patch
  object, the constructed faces will be added to the existing object.
  
`rStart`
: Coordinate(s) of the starting point, either a 3 element vector or
  a matrix with dimensions of $$[3\times n_{lineSegment}] to plot multiple line
  segments.
  
`rEnd`
: Coordinate(s) of the end point, either a 3 element vector or
  a matrix with dimensions of $$[3\times n_{lineSegment}]$$ to plot multiple line
  segments.
  
`r`
: Matrix with dimensions of $$[3\times n_{obj}\times n_{lineSegment}]$$. The function
  will plot $$n_{obj}$$ number of disconnected curves. The $$i$$th
  curve will follow the `x=r(1,i,:)`, `y=r(2,i,:)`, `z=r(3,i,:)`
  (parameteric) segmented curve.
  
`lineStyle`
: Line style, default value is `'-'` for continuous line. Any other
  Matlab line style string is accepted: `'--'`\|`':'`\|`'-.'`\|`'none'`.
  
`lineWidth`
: Line width in pt, default value is 0.5.
  
`mPatch`
: If `true`, a separate patch object will be created per line
  segment. Default is `false`, a single patch object will store all
  line segments.
  
### See Also
  
[line](https://www.mathworks.com/help/matlab/ref/line.html)
 

{% include links.html %}
