---
{title: swplot.figure, link: swplot.figure, summary: creates swplot figure, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_figure, folder: swplot, mathjax: 'true'}

---
  
### Syntax
  
`hFigure = swplot.figure`
 
`hFigure = swplot.figure(mode)`
  
### Description
  
`hFigure = swplot.figure` creates an empty figure with all the controls
for modifying the plot and the 3D roation engine initialized that rotates
the objects on the figure instead of the viewport. To plot anything onto
the figure, the handle of the graphics object (after creating it using
[surf](https://www.mathworks.com/help/matlab/ref/surf.html), [patch](https://www.mathworks.com/help/matlab/ref/patch.html), etc.) has to be added to the figure using
the function [swplot.add](swplot_add).
 
`hFigure = swplot.figure(mode)` defines settings for the figure.
 
### Input Arguments
  
`mode`
: Optional string. If `'nohg'`, then no hgtransform object will be
  used for fine object rotation. Can be usefull for certain
  export functions, that are incompatible with [hgtransform](https://www.mathworks.com/help/matlab/ref/hgtransform.html)
  objects. Default value is `'hg'` to use hgtransform.
  
### See Also
  
[swplot.add](swplot_add) \| [hgtransform](https://www.mathworks.com/help/matlab/ref/hgtransform.html)
 

{% include links.html %}
