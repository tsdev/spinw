---
{title: swplot.base( ), summary: sets the basis vectors of an swplot figure, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_base.html, folder: swplot, mathjax: 'true'}

---
sets the basis vectors of an [swplot figure](swplot_figure.html)
 
SWPLOT.BASE(BV, {hFigure})
 
BV is a matrix with dimensions of [3 3] and contains the three basis
vectors of the new coordinate system as column vectors.
 
SWPLOT.BASE(obj, {hFigure})
 
obj is a [spinw](spinw.html) object that defines the swplot coordinate system as
lattice units.
 
BV = SWPLOT.BASE
 
Returns the basis vectors stored in the [swplot figure](swplot_figure.html).
 
Input:
 
BV            Either a 3x3 matrix of the new basis vectors or a [spinw](spinw.html)
              object where the new basis vectors will be the lattice
              units of the stored crystal.
hFigure       Handle of the [swplot figure](swplot_figure.html). Default is the selected
              figure.
 
See also SWPLOT.PLOT.
 

