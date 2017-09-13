---
{title: swplot.legend( ), summary: draws legend to the swplot figure, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_legend.html, folder: swplot, mathjax: 'true'}

---
draws legend to the [swplot figure](swplot_figure.html)
 
SWPLOT.LEGEND({switch}, {hFigure})
 
status = SWPLOT.LEGEND
 
Input:
 
switch        One of the following string:
                  'on'                show legend,
                  'off'               hide legend,
                  'refresh'           redraw legend,
                  {'-','--','none'}   change the linestyle of the legend
                                      frame.
              Default is 'on'.
hFigure       Handle of the [swplot figure](swplot_figure.html). Default is the selected
              figure.
 
Example:
  [swplot.figure](swplot_figure.html)
  swplot.addcircle([0 0 0],[0 0 1],1)
  [swplot.legend](swplot_legend.html)
 

