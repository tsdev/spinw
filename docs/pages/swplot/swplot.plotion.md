---
{title: swplot.plotion, link: swplot.plotion, summary: plots magnetic ion properties,
  keywords: sample, sidebar: sw_sidebar, permalink: swplot_plotion.html, folder: swplot,
  mathjax: 'true'}

---

### Syntax

` `

### Description

 
The function plots selected properties of magnetic ions stored in a SpinW
object onto an swplot figure.
 

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
: String, defines how the bond is plotted
    'aniso'     Ellipsoid is plotted for single ion anisotropy.
    'g'     	Ellipsoid is drawn for g-tensor.

% `scale`
: Scaling factor for the size of the ellipsoid relative to the 
 hortest bond length. Default value is 1/3.

% `alpha`
: Transparency (alpha value) of the ellipsoid, default value is 
 .3.

% `radius1`
: Minimum radius of the ellipsoid, default value is 0.08.

% `lineWidth`
: Line width in pt of the main circles surrounding the ellipsoid, 
 f zero no circles are drawn. Default is 0.5.

% `figure`
: Handle of the swplot figure. Default is the selected figure.

% `legend`
: Whether to add the plot to the legend, default is true.

% `color`
: Color of the ellipsoid:
    'auto'      All ellipsoid gets the color of the ion.
    'colorname' All ellipsoid will have the same given color.
    [R G B]     RGB code of the color that fix the color of all
                ellipsoid.

% `color2`
: Color of the main circles, default is 'auto' when the ellipses
 ill have the same color as the ellipsoids. Can be either a row
 ector of RGB code or string of a color name.

% `nMesh`
: Resolution of the ellipse surface mesh. Integer number that is
 sed to generate an icosahedron mesh with #mesh number of
 dditional triangulation, default value is stored in
 wpref.getpref('nmesh')

% `nPatch`
: Number of points on the curve for the arrows, default
 alue is stored in swpref.getpref('npatch').

% `tooltip`
: If true, the tooltips will be shown when clicking on atoms.
 efault is true.

% `shift`
: Column vector with 3 elements, all atomic positions will be
 hifted by the given value. Default value is [0;0;0].

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
The name of the objects that are created called 'bond'. To find the
handles and the stored data on these objects, use e.g.
sObject = swplot.findobj(hFigure,'name','bond')

