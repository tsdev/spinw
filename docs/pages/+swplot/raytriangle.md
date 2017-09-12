---
{title: swplot.raytriangle( ), summary: finds if a ray crosses a triangle, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_raytriangle.html, folder: +swplot, mathjax: 'true'}

---
finds if a ray crosses a triangle
 
The code is optimised for a single ray.
 
SWPLOT.RAYTRIANGLE(V,F,ray)
 
Input:
 
V         Vertex positions in a matrix with dimensions [nVertex 3].
F         Faces in a matrix with dimensions [nFace 3], where 
              max(F) = nVertex.
ray       Definition of the ray via 2 points in space, while the ray
          pointing P1-->P2, stored in a matrix [P1;P2] with dimensions 
          [2 3].
 
