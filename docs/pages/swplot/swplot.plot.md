---
{title: swplot.plot, link: swplot.plot, summary: plots objects to swplot figure, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_plot, folder: swplot, mathjax: 'true'}

---
  
### Syntax
  
`swplot.plot(Name,Value)`
 
`hFigure = swplot.plot(Name,Value)`
  
### Description
  
`swplot.plot(Name,Value)` plots objects to the swplot figure and adds the
objects to the [hgtransform](https://www.mathworks.com/help/matlab/ref/hgtransform.html) object. This command enables the
plotting of multiple objects simultaneously while enabling fine control
of color, legend test, tooltip text etc. This commands is used by the
[spinw.plot](spinw_plot) high level plot command.
  
### Name-Value Pair Arguments
  
`'type'`
: Type of object to plot in a string. Possible options are:
  * `'arrow'`         position specifies start and end points,
  * `'ellipsoid'`     position specifies center,
  * `'cylinder'`      position specifies start and end points,
  * `'polyhedron'`    position specifies the vertices of the
                      convex polyhedron or polygon,
  * `'circle'`        position specifies center and normal vector,
  * `'line'`          position specifies start and end points (or
                      any number of points per curve),
  * `'text'`          position specifies the center of the text.
  
`'position'`
: Position of the object/objects in a matrix with dimensions of
  $$[3\times n_{obj}\times 2]$$ or $$[3\times n_{obj}]$$ or $$[3\times
  n_{obj}\times n_{point}]$$ depending on the type of object. The unit of
  the positions is determined by the `unit` parameter.
  
`'name'`
: String, the name of the object. It can be used for grouping the
  object handles to enable easier search, see [swplot.findobj](swplot_findobj) for
  details.
  
`'text'`
: Text to appear in the tooltip of the swplot figure after
  clicking on the object. Can be a string that will be the same
  for all objects, or a cell of strings for different text per
  object. Default value is taken from the label option.
  
`'label'`
: Text to appear in the legend in a string for the same text of
  all objects or strings in a cell with $$n_{obj}$$ number of elements for
  multiple objects. Default value is taken from the `name` parameter.
  
`'legend'`
: Type of legend to show the object:
  * `0`       do not show in legend,
  * `1`       colored box in legend,
  * `2`       dashed box in legend,
  * `3`       colored sphere in legend.
  
`'color'`
: Color of objects, either a single color or as many colors as
  many objects are given in a matrix with dimensions of $$[3\times 1]$$ or
  $$[3\times n_{obj}]$$. Colors are RGB triplets with values between 0 and
  255. Can be also string or cell of strings with the name of the colors,
  for possible color names see [swplot.color](swplot_color). Default value is `'red'`.
  
`'alpha'`
: Transparency of objects (1: non-transparent, 0: transparent)
  defined as a single number for uniform transparency or as a
  row vector with $$n_{obj}$$ number of elements to set transparency per object.
  Default value is 1.
  
`'unit'`
: String that determines the coordinate system where position vectors are
  defined:
  * `'lu'`    Lattice units are used where the lattice is defined
              by the stored basis (default).
  * `'xyz'`   Use the original Matlab units.
  
`'figure'`
: Handle of the swplot figure, default is the active figure.
  
`'R'`
: Radius value of cylinder, sphere (if no `'T'` parameter is given) and
  arrow, default value is 0.06.
  
`'ang'`
: Angle for arrow head in degree, default value is 15°.
  
`'lHead'`
: Length of the arrow head, default value is 0.5.
  
`'T'`
: Transformation matrix that transforms a unit sphere to the
  ellipse via: `R' = T(:,:,i)*R`, stored in a matrix with
  dimensions of $$[3\times 3\times n_{obj}]$$.
  
`'lineStyle'`
: Line style, default value is `'-'` for continuous lines. It can
  be also a vector with as many elements as many line segments.
  In this case the numbers are equivalent to the following style
  format string:
  * `1`   `'-'`,
  * `2`   `'--'`,
  * `3`   `'-.'`,
  * `4`   `':'`,
  * `5`   `'none'`.
  
`'lineWidth'`
: Line width, default value is 0.5, can be a vector with $$n_{obj}$$
  columns for different width per line segment.
  
`'fontSize'`
: Font size of text in pt when `type` parameter is set to `'text'`.
  Default value is stored in `swpref.getpref('fontsize')`.
  
`'nMesh'`
: Resolution of the ellipse surface mesh. Integer number that is
  used to generate an icosahedron mesh with `nMesh` number of
  additional subdivision of triangular surfaces. Default value is stored in
  `swpref.getpref('nmesh')`.
  
`'nPatch'`
: Number of points on the curve for arrow and cylinder, default
  value is stored in `swpref.getpref('npatch')`.
  
`'tooltip'`
: If `true`, the tooltip will be switched on after the
  plot. Default is `true`.
  
`'replace'`
: If `true`, all objects with the same name as the new plot will be
  deleted before plotting. Default is `false`.
  
`'data'`
: User supplied data per object that will be stored in the swplot
  figure and can be retrieved using [swplot.getdata](swplot_getdata). It is stored in a
  cell with $$n_{obj}$$ number of elements.
  
`'translate'`
: If `true`, the average center of the plot objects will be translated to
  the figure center. Default is `true`.
  
`'zoom'`
: If `true`, the swplot figure will be zoomed to make the plot objects
  cover the full figure. Default is `true`.
  
### See Also
  
[swplot.color](swplot_color) \| [swplot.add](swplot_add)
 

{% include links.html %}
