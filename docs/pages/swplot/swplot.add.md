---
{title: swplot.add( ), summary: adds a graphical object to the hgtransform of an swplot
    figure, keywords: sample, sidebar: sw_sidebar, permalink: swplot_add.html, folder: swplot,
  mathjax: 'true'}

---
adds a graphical object to the hgtransform of an [swplot figure](swplot_figure.html)
 
SWPLOT.ADD(hAdd, {hFigure})
 
It adds a graphical object to the hgtransform object of the figure to
enable continuous rotation with the mouse.
 
Input:
 
hAdd      Either vector of the handles of the graphical objects, or
          struct with dimensions of [1 nObject] with a handle field each
          contains a graphical object handle. The struct can contain any
          number of the following fields as well:
              'name'      Default value is 'general' if not given. The
                          name identifies groups of objects.
              'text'      Text that is shown in the tooltip when clicking
                          on the object.
              'position'  Position of the object, see [swplot.plot](swplot_plot.html) for
                          details.
              'label'     Label that is shown in the legend.
              'legend'    Type of legend, see [swplot.plot](swplot_plot.html) for details.
              'type'      Type of graphical object, see [swplot.plot](swplot_plot.html).
              'data'      Arbitrary data assigned to the object.
hFigure   The handle of the figure (number in the figure title bar). The
          default is the active [swplot figure](swplot_figure.html) if the argument is not
          provided by the user or it is empty matrix.
 
See also SWPLOT, SWPLOT.FIGURE, HGTRANSFORM, SWPLOT.PLOT.
 

