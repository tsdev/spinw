---
{title: swplot.icomesh, link: swplot.icomesh, summary: creates mesh by subdividing
    icosahedron faces, keywords: sample, sidebar: sw_sidebar, permalink: swplot_icomesh,
  folder: swplot, mathjax: 'true'}

---

### Syntax

`tr = swplot.icomesh(nsub)`

### Description

The output is a triangulated surface of the unit sphere, containing
20*4^nSub triangular faces. The output can be plotted using the trimesh()
function.
 

### Input Arguments

`nSub`
: Number of subdivisions. Default is 0 for icosahedron mesh
  output.

### Output Arguments

TR        Class triangulation object for plotting with trimesh().

{% include links.html %}
