---
{title: swplot.plotchem, link: swplot.plotchem, summary: plots polyhedra or chemical
    bonds, keywords: sample, sidebar: sw_sidebar, permalink: swplot_plotchem.html,
  folder: swplot, mathjax: 'true'}

---

### Syntax

` `

### Description

 
The function polyhedra around selected  atoms, or chemical bonds between
atoms an swplot figure.
 

### Input Arguments

### Name-Value Pair Arguments

% `obj`
: SpinW object.

% `range`
: Plotting range of the lattice parameters in lattice units,
 imensions are [3 2]. For example to plot the first unit cell,
 se: [0 1;0 1;0 1]. Also the number unit cells can be given
 long the a, b and c directions: [2 1 2], that is equivalent to
 0 2;0 1;0 2]. Default is the single unit cell.

% `unit`
: Unit in which the range is defined. It can be the following
 tring:
    'lu'        Lattice units (default).
    'xyz'       Cartesian coordinate system in Angstrom units.

% `mode`
: Selects the type of the plot:
    'poly'      Draws polyhedra around the center atoms
                (default).
    'bond'      Draws bonds between given atoms.

% `atom1`
: Indices of atoms stored in spinw.unit_cell for the center atom
 f the polyhedra or the first atom of the bonds. Can be also a
 tring that identifies the atoms by their labels.

% `atom2`
: Indices or label of the atoms stored in spinw.unit_cell. It
 etermines the vertices of the polyhedra or gives the second
 tom of a bond.

% `extend`
: If true, only atom1 has to be within the plotting range, atom2
 ill be searched +/-2 extra unit cells around the plotting
 ange. If false, both atom1 and atom2 have to be within the
 nit cell.

% `limit`
: Can be a single number which will restrict the number of
 earest negihbours of atom1 to connect. Can be also a vector
 hat defines bonds/polyhedra between atoms that are within the
 iven distance range stored as a row vector [dmin dmax].
 efault is 6 to plot octahedra around atom1.

% `alpha`
: Transparency of the plotted surfaces between 0 and 1 (1 for
 paque, 0 for transparent). Default value is 1 for bonds and
 .3 for polyhedron.

% `color`
: Surface color of the objects. Default is 'auto', when they are
 et to the color of atom1. [R G B] will fix the color of all
 onds to a uniform one, can be also arbitrary color name (see
 wplot.color() function). Can be also 'none', when no faces
 ill be shown.

% `color2`
: Color of the edges of the polyhedra (unused for bonds), default
 alue is 'auto' when the edge gets the same color as the faces.
 none' will remove the edges.

% `radius0`
: Radius of the cylinder, default value is 0.03.

% `figure`
: Handle of the swplot figure. Default is the selected figure.

% `legend`
: Whether to add the plot to the legend, default is true.

% `nPatch`
: Number of points on the curve for the cylinder, default
 alue is stored in swpref.getpref('npatch').

% `tooltip`
: If true, the tooltips will be shown when clicking on atoms.
 efault is true.

% `shift`
: Column vector with 3 elements, all atomic positions will be
 hifted by the given value by Angstrom units. Default value is
 0;0;0].

% `replace`
: Replace previous atom plot if true. Default is true.

% `translate`
: If true, all plot objects will be translated to the figure
 enter. Default is false.

% `zoom`
: If true, figure will be automatically zoomed to the ideal size.
 efault is false.

% `copy`
: If true, a hardcopy of the spinw object will be sved in the
 igure data, otherwise just the handle of the spinw object, 
 hus the figure can be updated when the spin object changed.
 efault value is false. 

### Output Arguments

hFigure           Handle of the swplot figure.
The name of the objects that are created called 'chem'. To find the
handles and the stored data on these objects, use e.g.
sObject = swplot.findobj(hFigure,'name','chem')

