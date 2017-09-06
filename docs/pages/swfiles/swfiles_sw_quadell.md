---
{title: sw_quadell( ), keywords: sample, summary: calculates and plots the parameters of an ellipsoid from a quadratic form,
  sidebar: sw_sidebar, permalink: swfiles_sw_quadell.html, folder: swfiles, mathjax: 'true'}

---
  calculates and plots the parameters of an ellipsoid from a quadratic form
 
  ellM = SW_QUADELL(M, plot)
 
  It calculates the parameters of an ellipsoid belonging to the input
  quadratic form and plots it. Can be used as a tool to visualise
  anisotropy matrices, g-tensors, etc.
 
  Input:
 
  M         Quadratic form coefficients in a 3x3 matrix.
  toplot    If true, the ellipsoid is plotted. Default is true.
 
  Output:
 
  ellM      Contais the principal axis of the ellipsoid in the columns, 3x3
            matrix.
 
  On the plot the principal axes are shown by a red line, the x,y and z
  axes are shown by black line.
 
