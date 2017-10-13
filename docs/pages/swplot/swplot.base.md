---
{title: swplot.base, link: swplot.base, summary: sets the basis vectors of an swplot
    figure, keywords: sample, sidebar: sw_sidebar, permalink: swplot_base, folder: swplot,
  mathjax: 'true'}

---
  
### Syntax
  
`swplot.base(BV)`
 
`swplot.base(obj)`
  
`BV = swplot.base`
 
### Description
  
`swplot.base(BV)` sets the basis vector for an swplot figure. The basis
vectors can be used to define a non-orthogonal coordinate system for
graphic objects.
 
`swplot.base(obj)` sets the basis vectors to the lattice units of a given
[spinw] object `obj`.
   
`BV = swplot.base` returns the basis vectors stored in the swplot figure.
   
  
### Input Arguments
  
`BV`
: Either a $$[3\times 3]$$ matrix of the new basis vectors or a [spinw]
  object where the new basis vectors will be the lattice
  units and the basis vectors are generated via [spinw.basisvector].
  
`hFigure`
: Handle of the [swplot] figure. Default is the active
  figure.
  
### See Also
  
[swplot.plot](swplot_plot)
 

{% include links.html %}
