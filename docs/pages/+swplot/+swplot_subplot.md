---
{title: subplot( ), keywords: sample, summary: create subplots with variable gaps between axes,
  sidebar: sw_sidebar, permalink: +swplot_subplot.html, folder: +swplot, mathjax: 'true'}

---
  create subplots with variable gaps between axes
 
  SWPLOT.SUBPLOT(m,n,p,space)
 
  SWPLOT.SUBPLOT([m n p],space)
 
  Input:
 
  m,n,p     Three integer that defines subplot, for details see the
            built-in subplot command.
  space     Vector with elements: [margin hgap vgap], where:
                margin  Top and right margin at the figure edge.
                hgap    Left margin and horizontal gap between axes.
                vgap    Bottom margin and vertical gap between axes.
 
  See also SUBPLOT.
 
