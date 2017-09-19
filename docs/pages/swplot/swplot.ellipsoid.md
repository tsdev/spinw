---
{title: swplot.ellipsoid, link: swplot.ellipsoid, summary: draw ellipsoid, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_ellipsoid.html, folder: swplot, mathjax: 'true'}

---

### Syntax

` `

### Description

Significant speedup can be achieved by a single patch compared to ellipse
per patch.
 
hPatch = SWPLOT.ELLIPSOID(handle,...
 
Handle can be the handle of an axes object or a patch object. It either
selects an axis to plot or a patch object (triangulated) to add vertices
and faces.
 

### Input Arguments

% `handle`
:   Handle of an axis or patch object. In case of patch object, the

% `constructed`
:cted faces will be added to the existing object instead

% `of`
:ting a new one.

% `R0`
:   Center of the ellipsoid stored in a matrix with dimensions of

% `[3`
:ipse].

% `T`
:   Transformation matrix that transforms a unit sphere to the

% `ellipse`
: via: R' = T(:,:,i)*R

% `Dimensions`
:ons are [3 3 nEllipse].

% `mesh`
:   Mesh of the ellipse surface, a triangulation class object or an

% `integer`
: that used to generate an icosahedron mesh with #mesh

% `number`
:of additional triangulation. Default value is stored in

% ``
:getpref('nmesh')

### See Also

[triangulation] and [swplot.icomesh](swplot_icomesh.html)

