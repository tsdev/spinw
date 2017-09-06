---
{title: polyhedron( ), keywords: sample, summary: draw convex polyhedra or polynom from vertex list,
  sidebar: sw_sidebar, permalink: +swplot_polyhedron.html, folder: +swplot, mathjax: 'true'}

---
  draw convex polyhedra or polynom from vertex list
 
  hPatch = SWPLOT.POLYHEDRON(vertices)
 
  hPatch = SWPLOT.POLYHEDRON(handle,...
 
  Handle can be the handle of an axes object or a patch object. It either
  selects an axis to plot or a patch object (triangulated) to add vertices
  and faces.
 
  Input:
 
  vertices      Matrix with dimensions [3 nObject nPoint], where nObject is
                the number of polyhedra to draw, nPoint is the number of
                vertices per polyhedron.
 
  See also CONVHULLN.
 
