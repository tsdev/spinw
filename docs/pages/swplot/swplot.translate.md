---
{title: swplot.translate, link: swplot.translate, summary: translate objects on swplot
    figure, keywords: sample, sidebar: sw_sidebar, permalink: swplot_translate, folder: swplot,
  mathjax: 'true'}

---

### Syntax

`swplot.translate(mode, {hfigure})`

### Description

The function translates the objects of an swplot, where the coordinate
system is defined by the plane of the figure with horizontal x-axis,
vertical y-axis and out-of-plane z-axis.
 

### Input Arguments

`mode`
: Either a vector with three elements determining the translation 
  value in the figure plane coordinate system, or 'auto' that
  centers the object on the figure. Default is 'auto'.

`hFigure`
: Handle of the swplot figure window, optional.

### See Also

[swplot.zoom](swplot_zoom)

{% include links.html %}
