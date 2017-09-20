---
{title: swplot.polyhedron, link: swplot.polyhedron, summary: draw convex polyhedra
    or polynom from vertex list, keywords: sample, sidebar: sw_sidebar, permalink: swplot_polyhedron.html,
  folder: swplot, mathjax: 'true'}

---

### Syntax

`hpatch = swplot.polyhedron(vertices)`

### Description

hPatch = SWPLOT.POLYHEDRON(handle,...
 
Handle can be the handle of an axes object or a patch object. It either
selects an axis to plot or a patch object (triangulated) to add vertices
and faces.
 

### Input Arguments

`vertices`
: Matrix with dimensions [3 nObject nPoint], where nObject is
  the number of polyhedra to draw, nPoint is the number of
  vertices per polyhedron.

### See Also

[convhulln]

