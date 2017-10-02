---
{title: swplot.plotmag, link: swplot.plotmag, summary: plots magnetic structure, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_plotmag, folder: swplot, mathjax: 'true'}

---

### Syntax

`swplot.plotmag(Name,Value)`

### Description

hFigure = SWPLOT.PLOTMAG(...)
 
The function plots the magnetic structure of a SpinW object onto an
swplot figure.
 

### Input Arguments

### Name-Value Pair Arguments

`'obj'`
: SpinW object.

`'range'`
: Plotting range of the lattice parameters in lattice units,
  dimensions are [3 2]. For example to plot the first unit cell,
  use: [0 1;0 1;0 1]. Also the number unit cells can be given
  along the a, b and c directions: [2 1 2], that is equivalent to
  [0 2;0 1;0 2]. Default is the single unit cell.

`'unit'`
: Unit in which the range is defined. It can be the following
  string:
      'lu'        Lattice units (default).
      'xyz'       Cartesian coordinate system in Å units.

`'mode'`
: String, defines the way the magnetic moments are plotted:
      'all'       Plot both the rotation plane of incommensurate
                  magnetic structures and the moment directions.
      'circle'    Plot only the rotation plane of incommensurate
                  magnetic structures.
      'arrow'     Plots only the moment directions.
      'none'      Don't plot anything.

`'figure'`
: Handle of the swplot figure. Default is the selected figure.

`'legend'`
: Whether to add the plot to the legend, default value is true.

`'label'`
: Whether to plot labels for atoms, default value is true.

`'dText'`
: Distance between item and its text label, default value is 0.1
  Å.

`'fontSize'`
: Font size of the atom labels in pt, default value is stored in
  swpref.getpref('fontsize').

`'color'`
: Color of the magnetic moments:
      'auto'      All moments get the same color as the magnetic
                  atom.
      'colorname' All moments will have the same color.
      [R G B]     RGB code of the color.

`'scale'`
: Scaling factor for the lenght of the magnetic moments relative
  to the length of the shortest bond (if there are no bonds, 3A 
  is taken as bond length). Default is 0.4.

`'normalize'`
: If true, all moment length will be normalized to the scale
  factor, default value is false.

`'radius0'`
: Radius value of arrow body, default value is 0.06.

`'ang'`
: Angle of the arrow head in ° units, default value is 30 °.

`'lHead'`
: Length of the arrow head, default value is 0.5.

`'α'`
:   Transparency (α value) of the circle, representing the
  rotation plane of the moments, default value is 0.07.

`'centered'`
: If true, the moment vector is centered on the atom, if false
  the beggining of the spin vector is on the atom. Default is
  true.

`'nPatch'`
: Number of points on the curve for the arrows, default
  value is stored in swpref.getpref('npatch').

`'tooltip'`
: If true, the tooltips will be shown when clicking on atoms.
  Default is true.

`'shift'`
: Column vector with 3 elements, all vectors will be
  shifted by the given value. Default value is [0;0;0].

`'replace'`
: Replace previous magnetic moment plot if true. Default is true.

`'translate'`
: If true, all plot objects will be translated to the figure
  center. Default is false.

`'zoom'`
: If true, figure will be automatically zoomed to the ideal size.
  Default is false.

`'copy'`
: If true, a hardcopy of the spinw object will be sved in the
  figure data, otherwise just the handle of the spinw object, 
  thus the figure can be updated when the spin object changed.
  Default value is false. 

{% include links.html %}
