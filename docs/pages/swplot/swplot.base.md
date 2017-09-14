---
{title: swplot.base, link: swplot.base, summary: sets the basis vectors of an swplot
    figure, keywords: sample, sidebar: sw_sidebar, permalink: swplot_base.html, folder: swplot,
  mathjax: 'true'}

---
 
SWPLOT.BASE(BV, {hFigure})
 
BV is a matrix with dimensions of [3 3] and contains the three basis
vectors of the new coordinate system as column vectors.
 
SWPLOT.BASE(obj, {hFigure})
 
obj is a spinw object that defines the swplot coordinate system as
lattice units.
 
BV = SWPLOT.BASE
 
Returns the basis vectors stored in the swplot figure.
 
Input:
 
BV            Either a 3x3 matrix of the new basis vectors or a spinw
              object where the new basis vectors will be the lattice
              units of the stored crystal.
hFigure       Handle of the swplot figure. Default is the selected
              figure.
 
See also SWPLOT.PLOT.
 

