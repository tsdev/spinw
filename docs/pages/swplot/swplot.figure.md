---
{title: swplot.figure, link: swplot.figure, summary: creates swplot figure, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_figure.html, folder: swplot, mathjax: 'true'}

---
 
hFigure = SWPLOT.FIGURE({mode})
 
The function creates an empty figure with all the controls for modifying
the plot and the 3D roation engine that rotates the objects on the figure
instead of the viewport. To plot anything onto the figure, the handle of
the graphics object (after creating it using surf, patch, etc.) has to be
added to the list using the function swplot.add.
 
Input:
 
mode      Optional string. If 'nohg', then no hgtransform object will be
          used for fine object rotation. Can be usefull for certain
          export functions, that are incompatible with hgtransform
          objects. Default is 'hg' to use hgtransform.
 
See also SWPLOT.ADD, HGTRANSFORM.
 

