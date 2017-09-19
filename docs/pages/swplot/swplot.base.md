---
{title: swplot.base, link: swplot.base, summary: sets the basis vectors of an swplot
    figure, keywords: sample, sidebar: sw_sidebar, permalink: swplot_base.html, folder: swplot,
  mathjax: 'true'}

---

### Syntax

` `

### Description

vectors of the new coordinate system as column vectors.
 
SWPLOT.BASE(obj, {hFigure})
 
obj is a spinw object that defines the swplot coordinate system as
lattice units.
 
BV = SWPLOT.BASE
 
Returns the basis vectors stored in the swplot figure.
 

### Input Arguments

% `BV`
:     Either a 3x3 matrix of the new basis vectors or a spinw

% `object`
:ere the new basis vectors will be the lattice

% `units`
:the stored crystal.

% `hFigure`
:     Handle of the swplot figure. Default is the selected

% ``
:

### See Also

[swplot.plot](swplot_plot.html)

