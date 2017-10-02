---
{title: swplot.plotatom, link: swplot.plotatom, summary: plots crystal structure,
  keywords: sample, sidebar: sw_sidebar, permalink: swplot_plotatom, folder: swplot,
  mathjax: 'true'}

---

### Syntax

`swplot.plotatom(Name,Value)`

### Description

hFigure = SWPLOT.PLOTATOM(...)
 
The function plots the crystal structure of a SpinW object onto an swplot
figure.
 

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
: String, defines the types of atoms to plot:
      'all'       Plot all atoms (default).
      'mag'       Plot magnetic atoms only.
      'nonmag'    Plot non-magnetic atoms only.

`'figure'`
: Handle of the swplot figure. Default is the selected figure.

`'legend'`
: Whether to add the plot to the legend, default value is true.

`'label'`
: Whether to plot labels for atoms, default value is false.

`'dText'`
: Distance between item and its text label, default value is 0.1
  Å.

`'fontSize'`
: Font size of the atom labels in pt, default value is stored in
  swpref.getpref('fontsize').

`'radius0'`
: Constant atom radius, default value is 0.3 Å.

`'radius'`
: Defines the atom radius:
      'fix'       Sets the radius of all atoms to the value
                  stored in radius0.
      'auto'      use radius data from database based on the atom
                  label multiplied by radius0 value.

`'color'`
: Color of the atoms:
      'auto'      All atom gets the color stored in obj.unit_cell.
      'colorname' All atoms will have the same color.
      [R G B]     RGB code of the color that fix the color of all
                  atoms.

`'nMesh'`
: Resolution of the ellipse surface mesh. Integer number that is
  used to generate an icosahedron mesh with #mesh number of
  additional triangulation, default value is stored in
  swpref.getpref('nmesh')

`'tooltip'`
: If true, the tooltips will be shown when clicking on atoms.
  Default is true.

`'shift'`
: Column vector with 3 elements, all atomic positions will be
  shifted by the given value in Å units. Default value is
  [0;0;0].

`'replace'`
: Replace previous atom plot if true. Default is true.

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

### Output Arguments

hFigure           Handle of the swplot figure.
The name of the objects that are created called 'atom' and 'atom_label'.
To find the handles and the stored data on these objects, use e.g.
  sObject = swplot.findobj(hFigure,'name','atom')

{% include links.html %}
