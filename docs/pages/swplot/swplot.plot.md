---
{title: swplot.plot, link: swplot.plot, summary: plots objects to swplot figure, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_plot.html, folder: swplot, mathjax: 'true'}

---

### Syntax

` `

### Description

 
 

### Name-Value Pair Arguments

% `type`
:   Type of object to plot in a string. Possible options are:
 row'         position specifies start and end points
 lipsoid'     position specifies center
 linder'      position specifies start and end points
 lyhedron'    position specifies the vertices of the
              convex polyhedron or polygon
 rcle'        position specifies center and normal vector
 ne'          position specifies start and end points (or
              any number of points per curve)
 xt'          position specifies the center of the text

% `position`
:n  Position of the object/objects in a matrix with dimensions of

% `[3`
:ect 2]/[3 nObject]/[3 nObject nPoint] depending on the

% `type`
: object. The unit of the positions is determined by the

% `'unit'`
:option.

% `name`
:   String, the name of the object. It can be used for finding the

% `object`
:handles after plotting.

% `text`
:   Text to appear in the tooltip of the swplot figure after

% `clicking`
:g on the object. Can be a string that will be the same

% `for`
: objects, or a cell of strings for different text per

% `object.`
: Default value is taken from the label option.

% `label`
:   Text to appear in the legend in a string for the same text of

% `all`
:ects or strings in a cell for multiple objects with

% `dimension`
:on [1 nObject]. Default value is taken from the name

% ``
:

% `legend`
:   Type of legend to show the object:
      Do not show in legend.
      Colored box in legend.
      Dashed box in legend.
      Colored sphere in legend.

% `color`
:   Color of objects, either a single color or as many colors as

% `many`
:jects are given in a matrix with dimensions of [3 1]/[3

% `nObject].`
:]. Values are RGB triples with values between [0 255].

% `Can`
:also string or cell of strings with the name of the

% `colors,`
: for details see swplot.color. Default is red.

% `alpha`
:   Transparency of objects (1 non-transparent, 0 transparent)

% `defined`
: as a single number for unitform transparency or as a

% `row`
:tor with nObject element to set transparency per object.

% `Default`
: value is 1.

% `unit`
:   String that determines the coordinate system:
 '    Lattice units are used where the lattice is defined
      by the stored basis (default).
 z'   Use the original matlab units.

% `figure`
:   Handle of the swplot figure. Default is the selected figure.

% `R`
:   Radius value of cylinder, sphere (if no 'T' is given) and

% `arrow,`
:default is 0.06.

% `ang`
:   Angle for arrow head in degree units, default is 15 degree.

% `lHead`
:   Length of the arrow head, default value is 0.5.

% `T`
:   Transformation matrix that transforms a unit sphere to the

% `ellipse`
: via: R' = T(:,:,i)*R

% `Dimensions`
:ons are [3 3 nObject].

% `lineStyle`
:le Line style, default value is '-' for continuous lines. It can

% `be`
: a vector with as many elements as many line segments.

% `In`
: case the numbers are equivalent to the following style

% `format`
:string:
  '-'
  '--'
  '-.'
  ':'
  'none'

% `lineWidth`
:th Line width, default value is 0.5, can be a vector with nObject

% `columns`
: for different width per line segment.

% `fontSize`
:e  Font size of text when type option is set to 'text'. Default

% `value`
:s stored in swpref.getpref('fontsize').

% `nMesh`
:   Resolution of the ellipse surface mesh. Integer number that is

% `used`
: generate an icosahedron mesh with #mesh number of

% `additional`
:nal triangulation, default value is stored in

% ``
:getpref('nmesh')

% `nPatch`
:   Number of points on the curve for arrow and cylinder, default

% `value`
:s stored in swpref.getpref('npatch').

% `tooltip`
:   If true, the tooltip will be switched on at the end of the

% `plot.`
:efault is true.

% `replace`
:   If true, all object with the same name as the new plot will be

% `deleted`
: before plotting. Default is false.

% `data`
:   Arbitrary data per object that will be stored in the swplot

% `figure`
:and can be retrieved. It is stored in a cell with

% `nObject`
: number of elements.

% `translate`
:te If true, all plot objects will be translated to the figure

% `center.`
: Default is true.

% `zoom`
:   If true, figure will be automatically zoomed to the ideal size.

% `Default`
: is true.

### See Also

[swplot.color](swplot_color.html)

