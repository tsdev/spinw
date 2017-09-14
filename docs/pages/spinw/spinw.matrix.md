---
{title: spinw.matrix property, link: spinw.matrix, summary: Stores 3x3 matrices for
    using them in the Hailtonian., keywords: sample, sidebar: sw_sidebar, permalink: spinw_matrix.html,
  folder: spinw, mathjax: 'true'}

---
Sub fields are:
  mat     stores the actual values of 3x3 matrices, in a
          3 x 3 x nMatrix matrix, defult unit is meV
  color   color assigned for every matrix, stored in a
          3 x nMatrix matrix, with 0-255 RGB columns
  label   label for every matrix, stored as string in a
          1 x nMatrix cell
 
See also SPINW.ADDMATRIX, SPINW.NTWIN.

