---
{title: swplot.subfigure, link: swplot.subfigure, summary: changes position of figure
    window on the screen, keywords: sample, sidebar: sw_sidebar, permalink: swplot_subfigure.html,
  folder: swplot, mathjax: 'true'}

---
 
SWPLOT.SUBFIGURE(m,n,p,{hFigure})
 
Input:
 
Changes the position of the figure window on the screen, the position is
determined similarly to the Matlab function subplot(). Here the screen is
the canvas where the figure window is positioned.
 
SUBFIGURE divides the display into an m-by-n grid and moves the figure
window in the position specified by p. It numbers its figures by row,
such that the first figure is the first column of the first row, the
second figure is the second column of the first row, and so on.
 
Input:
 
m,n,p     Integer numbers, defining figure window position.
hFigure   Handle of the figure window, optional. Default is gcf.
 

