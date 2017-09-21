---
{title: sw_quadell, link: sw_quadell, summary: calculates and plots the parameters
    of an ellipsoid from a quadratic form, keywords: sample, sidebar: sw_sidebar,
  permalink: sw_quadell.html, folder: swfiles, mathjax: 'true'}

---

### Syntax

`ellm = sw_quadell(m, plot)`

### Description

It calculates the parameters of an ellipsoid belonging to the input
quadratic form and plots it. Can be used as a tool to visualise
anisotropy matrices, g-tensors, etc.
 

### Input Arguments

`M`
: Quadratic form coefficients in a 3x3 matrix.

`toplot`
: If true, the ellipsoid is plotted. Default is true.

### Output Arguments

ellM      Contais the principal axis of the ellipsoid in the columns, 3x3
          matrix.
On the plot the principal axes are shown by a red line, the x,y and z
axes are shown by black line.

{% include links.html %}
