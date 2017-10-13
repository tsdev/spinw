---
{title: swplot.ellipsoid, link: swplot.ellipsoid, summary: creates a 3D ellipsoid
    patch, keywords: sample, sidebar: sw_sidebar, permalink: swplot_ellipsoid, folder: swplot,
  mathjax: 'true'}

---
  
### Syntax
  
`hPatch = swplot.ellipsoid(R0, T)`
 
`hPatch = swplot.ellipsoid(R0, T, nMesh)`
  
`hPatch = swplot.ellipsoid(handle, ...)`
 
### Description
  
`hPatch = swplot.ellipsoid(R0,T)` creates multiple ellipsoids with a
single [patch](https://www.mathworks.com/help/matlab/ref/patch.html) command. The ellipsoids are defined by the position
of the center and a $$[3\times 3]$$ matrix, a qudratic form.
 
Significant speedup can be achieved using a single patch command to
generate many ellipsoids compared to drawing single ellipse per patch.
   
`hPatch = swplot.ellipsoid(R0, T, nMesh)` defines the size of the mesh
that defines the surface.
 
`hPatch = swplot.ellipsoid(handle, ...)` adds the generated patch object
to a given axis if `handle` is an axis handle or adds the ellipsoids to
an existing [patch](https://www.mathworks.com/help/matlab/ref/patch.html) object, if the given `handle` points to a
patch object.
  
### Input Arguments
  
`handle`
: Handle of an axis or triangulated patch object. In case of patch
  object, the constructed faces will be added to the existing object.
  
`R0`
: Center of the ellipsoids stored in a column vector with 3 elements or a
  matrix with dimensions of $$[3\times n_{obj}]$$ when multiple ellipsoids
  are defined at once.
  
`T`
: Transformation matrix that transforms a unit sphere to the desired
  ellipsoid by applying: `R' = T(:,:,i)*R`. In case of multiple
  ellipsoids the parameter is stored in a matrix with dimensions of
  $$[3\times 3\times n_{obj}]$$.
  
`nMesh`
: Mesh of the ellipse surface, a triangulation class object or an
  integer that used to generate an icosahedron mesh with $$n_{mesh}$$
  number of additional subdivision into triangles. Default value is stored in
  `swpref.getpref('nmesh')`, see also [swplot.icomesh](swplot_icomesh).
  
### See Also
  
[triangulation](https://www.mathworks.com/help/matlab/ref/triangulation.html) \| [swplot.icomesh](swplot_icomesh)
 

{% include links.html %}
