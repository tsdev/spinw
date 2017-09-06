---
{title: plot( ), keywords: sample, summary: 'plots crystal structure, magnetic structure,
    anisotropy and bonds', sidebar: sw_sidebar, permalink: '@spinw_plot.html', folder: '@spinw',
  mathjax: 'true'}

---
  plots crystal structure, magnetic structure, anisotropy and bonds
 
  PLOT(obj, 'Option1', Value1, ...)
 
  hFigure = PLOT(obj, 'Option1', Value1, ...)
 
  The function plots the atoms and couplings stored in obj onto an swplot
  figure. The figure can be rotated in 3D using the mouse and panning works
  by keeping the Ctrl/Control button pressed. There is information about
  every object on the figure (called tooltip) that is shown when clicked on
  the object. The 3D view direction can be changed programatically using
  the swplot.view() command while translations are controlled using the
  swplot.translate() command. Arbitrary transformation (combination of
  rotation and translation) can be introduced using the swplot.transform()
  command.
 
  The function calls several plot routine to draw different group of
  objects to the figure: atom (atoms), mag (magnetic moments), ion (single
  ion properties), bond (bonds), base (basis vectors) and cell (unit
  cells).
 
  Each group has a corresponding plot function with its own input options,
  for the allowed options, please check the help of the functions below.
 
  atom  swplot.plotatom
  mag   swplot.plotmag
  ion   swplot.plotion
  bond  swpplot.plotbond
  base  swplot.plotbase
  cell  swplot.plotcell
 
  To provide an option to any of these sub functions add the name of the
  group to the option. For example to set the 'color' option of the cell
  (change the color of the unit cell) use the option 'cellColor'.
 
  It is possible to switch of any of the subfunctions by using the option
  [groupName 'mode'] set to 'none'. For example to skip plotting the atoms
  use the 'atomMode' option set to 'none'.
 
  There are also global options, that each group recognizes, these global
  options can be added without the group name.
 
  Global options:
 
  range     Plotting range of the lattice parameters in lattice units,
            dimensions are [3 2]. For example to plot the first unit cell,
            use: [0 1;0 1;0 1]. Also the number unit cells can be given
            along the a, b and c directions: [2 1 2], that is equivalent to
            [0 2;0 1;0 2]. Default is the single unit cell.
  unit      Unit in which the range is defined. It can be the following
            string:
                'lu'        Lattice units (default).
                'xyz'       Cartesian coordinate system in Angstrom units.
  figure    Handle of the swplot figure. Default is the selected figure.
  legend    Whether to add the plot to the legend, default is true.
  fontSize  Font size of the atom labels in pt, default is stored in
            swpref.getpref('fontsize').
  nMesh     Resolution of the ellipse surface mesh. Integer number that is
            used to generate an icosahedron mesh with #mesh number of
            additional triangulation, default value is stored in
            swpref.getpref('nmesh')
  nPatch    Number of points on the curve for the arrows and cylinders,
            default value is stored in swpref.getpref('npatch').
  tooltip   If true, the tooltips will be shown when clicking on the plot.
            Default is true.
  shift     Column vector with 3 elements, all objects will be shifted by
            the given value. Default value is [0;0;0].
  replace   Replace previous plots if true. Default is true.
  translate If true, all plot objects will be translated to the figure
            center. Default is true.
  zoom      If true, figure will be automatically zoomed to the ideal size.
            Default is true.
 
