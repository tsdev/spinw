---
{title: swplot.plotbase( ), summary: plots the edges of unit cells on swplot figure,
  keywords: sample, sidebar: sw_sidebar, permalink: swplot_plotbase.html, folder: swplot,
  mathjax: 'true'}

---
plots the edges of unit cells on swplot figure
 
SWPLOT.PLOTBASE('Option1',Value1,...)
 
hFigure = SWPLOT.PLOTBASE('Option1',Value1,...)
 
Options:
 
mode      String that determines the type base to plot. Possible values
          are:
              abc     Plots the lattice vectors, default.
              xyz     Plots the lattice Descartes coordinate system.
              hkl     Plots the reciprocal lattice vectors.
length    Determines the length of the a, b and c arrows. If 0, the
          length will be equal to the corresponding lattice parameters,
          while if non-zero, the number determines the length in
          Angstrom. Default is 2 Angstrom.
label     Logical variable, plots abc labels if true. Default is true.
figure    Handle of the swplot figure. Default is the selected figure.
color     Color of the arrows, default is red-green-blue for abc, stored
          in the columns of a 3x3 matrix.
R         Radius value of arrow body, default is 0.06.
alpha     Head angle of the arrow in degree units, default is 30 degree.
lHead     Length of the arrow head, default value is 0.5.
d         Distance from origin in xyz units.
dtext     Distance from arrow and in xyz units.
shift     Column vector with 3 elements, the basis vectors will be
          shifted by the given values in Angstrom unit. Default value is
          [0;0;0].
translate If true, all plot objects will be translated to the figure
          center. Default is false.
zoom      If true, figure will be automatically zoomed to the ideal size.
          Default is false.
 
