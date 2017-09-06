---
{title: plotion( ), keywords: sample, summary: plots magnetic ion properties, sidebar: sw_sidebar,
  permalink: +swplot_plotion.html, folder: +swplot, mathjax: 'true'}

---
  plots magnetic ion properties
 
  SWPLOT.PLOTION('option1', value1, ...)
 
  hFigure = SWPLOT.PLOTION(...)
 
  The function plots selected properties of magnetic ions stored in a SpinW
  object onto an swplot figure.
 
  Input:
 
  Options:
 
  obj       SpinW object.
  range     Plotting range of the lattice parameters in lattice units,
            dimensions are [3 2]. For example to plot the first unit cell,
            use: [0 1;0 1;0 1]. Also the number unit cells can be given
            along the a, b and c directions: [2 1 2], that is equivalent to
            [0 2;0 1;0 2]. Default is the single unit cell.
  unit      Unit in which the range is defined. It can be the following
            string:
                'lu'        Lattice units (default).
                'xyz'       Cartesian coordinate system in Angstrom units.
  mode      String, defines how the bond is plotted
                'aniso'     Ellipsoid is plotted for single ion anisotropy.
                'g'     	Ellipsoid is drawn for g-tensor.
  scale     Scaling factor for the size of the ellipsoid relative to the 
            shortest bond length. Default value is 1/3.
  alpha     Transparency (alpha value) of the ellipsoid, default value is 
            0.3.
  radius1   Minimum radius of the ellipsoid, default value is 0.08.
  lineWidth Line width in pt of the main circles surrounding the ellipsoid, 
            if zero no circles are drawn. Default is 0.5.
  figure    Handle of the swplot figure. Default is the selected figure.
  legend    Whether to add the plot to the legend, default is true.
  color     Color of the ellipsoid:
                'auto'      All ellipsoid gets the color of the ion.
                'colorname' All ellipsoid will have the same given color.
                [R G B]     RGB code of the color that fix the color of all
                            ellipsoid.
  color2    Color of the main circles, default is 'auto' when the ellipses
            will have the same color as the ellipsoids. Can be either a row
            vector of RGB code or string of a color name.
  nMesh     Resolution of the ellipse surface mesh. Integer number that is
            used to generate an icosahedron mesh with #mesh number of
            additional triangulation, default value is stored in
            swpref.getpref('nmesh')
  nPatch    Number of points on the curve for the arrows, default
            value is stored in swpref.getpref('npatch').
  tooltip   If true, the tooltips will be shown when clicking on atoms.
            Default is true.
  shift     Column vector with 3 elements, all atomic positions will be
            shifted by the given value. Default value is [0;0;0].
  replace   Replace previous atom plot if true. Default is true.
  translate If true, all plot objects will be translated to the figure
            center. Default is false.
  zoom      If true, figure will be automatically zoomed to the ideal size.
            Default is false.
  copy      If true, a hardcopy of the spinw object will be sved in the
            figure data, otherwise just the handle of the spinw object, 
            thus the figure can be updated when the spin object changed.
            Default value is false. 
 
  Output:
 
  hFigure           Handle of the swplot figure.
 
  The name of the objects that are created called 'bond'. To find the
  handles and the stored data on these objects, use e.g.
 
    sObject = swplot.findobj(hFigure,'name','bond')
 
