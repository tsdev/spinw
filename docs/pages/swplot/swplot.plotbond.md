---
{title: swplot.plotbond, link: swplot.plotbond, summary: plots magnetic bonds, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_plotbond.html, folder: swplot, mathjax: 'true'}

---

### Syntax

` `

### Description

 
The function plots the magnetic bonds stored in a SpinW object onto an
swplot figure.
 

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
    'cylinder'  Bonds are plotted as cylinders (default).
    'arrow'     Bonds are plotted as arrows (default if DM
                interactions are non-zero).
    'line'      Bonds are plotted as lines.
    'empty'     No bonds are plotted.

% `mode2`
: String, defines what is plotted on the bond:
    'none'      Don't plot anything on the bond (default).
    'antisym'   Plot the antisymmetric part (DM vector) of the 
                exchange at the middle point of the bond
                (default if DM vectors are non-zero).
    'sym'       Plot the symmetric exchange at the middle
                of the bond as an ellipsoid.

% `sign`
: String, defines how the ellipsoids are generated for exchange
 atrices that contain both negative and positive eigenvalues.
 ossible values are:
    'abs'       The absolute value of the eigenvalues is used.
                This works nicely except that AFM and FM values
                will have the same radius. Default value.
    'min'       If there is a negative eigenvalue, it is
                shifted to zero with all other egeinvalues
                equally. This works nicely to emphasize AFM
                type values in the exchange matrix. Problem is
                that 0 exchange values can be assigned non-zero
                radius.
    'max'       Same as min, just positive eigenvalues are
                shifted to zero. This works nicely to emphasize
                FM type exchange values, has the same problem
                as the 'min' option.

% `linewidth`
: Defines the bond radius if it is drawn by a line:
    'fix'       All line will have a width given by linewidth0.
                Default value.
    'lin'       Lines will have a width that is depending 
                linearly on the exchange matrix on the bond:
                        Width ~ sum(abs(J)), 
                where the largest line width on
                the strongest bond is given by linewidth0.
    'pow'       Same as 'auto', but the line width is a
                power function of J: W~(sum(abs(J))).^widthpow

% `widthpow`
: Defines the power that determines the linewidth if 'linewidth'
 ption is 'pow'.

% `linewidth0`
:0 Line width in pt used to draw the bond if 'mode' is 'line'. 
 efault value is 0.5.

% `lineStyle`
: Determines the line style when bonds plotted as lines. Possible
 alues:
    'auto'      Bonds are plotted as continuous/dashed lines
                depending on the label of the corresponding
                matrix (dashed line is used if the matrix
                label ends with '-', otherwise continuous).
                Default value.
    '--'        Bonds are plotted as dashed lines.
    '-'         Bonds are plotted as lines.

% `zero`
: If true, bonds with zero exchange matrix will be plotted as
 ell. Default is true.

% `radius0`
: Radius of the cylinder, default value is 0.05.

% `radius1`
: Radius of the DM vector and the minimum radius of the 
 llipsoid, default value is 0.08.

% `radius2`
: Constant atom radius, default value is 0.3 Angstrom.

% `radius`
: Defines the atom radius (important for arrow bonds, to avoid
 verlap with the spheres of the atoms):
    'fix'       Sets the radius of all atoms to the value
                given by radius2.
    'auto'      use radius data from database based on the atom
                label multiplied by radius2 value.

% `ang`
: Angle of the arrow head in degree units, default is 30 degree.

% `lHead`
: Length of the arrow head, default value is 0.3.

% `scale`
: Scaling factor for the length of the DM vector or the size of
 he ellipsoid relative to the shortest bond length. Default 
 alue is 1/3.

% `figure`
: Handle of the swplot figure. Default is the selected figure.

% `legend`
: Whether to add the plot to the legend, default is true.

% `color`
: Color of the bonds:
    'auto'      All bonds get the stored color.
    'colorname' All bonds will have the same given color.
    [R G B]     RGB code of the color that fix the color of all
                bonds.

% `color2`
: Color of the ellipse or DM vector on the bond:
    'auto'      All object get the color of the bond.
    'colorname' All object will have the same given color.
    [R G B]     RGB code of the color that fix the color of all
                object.

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

