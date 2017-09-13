---
{title: swplot.icomesh( ), summary: creates mesh by subdividing icosahedron faces,
  keywords: sample, sidebar: sw_sidebar, permalink: swplot_icomesh.html, folder: swplot,
  mathjax: 'true'}

---
creates mesh by subdividing icosahedron faces
 
TR = SWPLOT.ICOMESH(nSub)
 
The output is a triangulated surface of the unit sphere, containing
20*4^nSub triangular faces. The output can be plotted using the trimesh()
function.
 
Input:
 
nSub      Number of subdivisions. Default is 0 for icosahedron mesh
          output.
 
Output:
 
TR        Class triangulation object for plotting with trimesh().
 

