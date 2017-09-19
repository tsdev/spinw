---
{title: swplot.add, link: swplot.add, summary: adds a graphical object to the hgtransform
    of an swplot figure, keywords: sample, sidebar: sw_sidebar, permalink: swplot_add.html,
  folder: swplot, mathjax: 'true'}

---

### Syntax

` `

### Description

enable continuous rotation with the mouse.
 

### Input Arguments

% `hAdd`
:  Either vector of the handles of the graphical objects, or

% `struct`
:ith dimensions of [1 nObject] with a handle field each

% `contains`
: a graphical object handle. The struct can contain any

% `number`
:f the following fields as well:
 e'      Default value is 'general' if not given. The
         name identifies groups of objects.
 t'      Text that is shown in the tooltip when clicking
         on the object.
 ition'  Position of the object, see swplot.plot for
         details.
 el'     Label that is shown in the legend.
 end'    Type of legend, see swplot.plot for details.
 e'      Type of graphical object, see swplot.plot.
 a'      Arbitrary data assigned to the object.

% `hFigure`
:  The handle of the figure (number in the figure title bar). The

% `default`
:is the active swplot figure if the argument is not

% `provided`
: by the user or it is empty matrix.

### See Also

[swplot], [swplot.figure](swplot_figure.html), [hgtransform] and [swplot.plot](swplot_plot.html)

