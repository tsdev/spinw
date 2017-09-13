---
{title: spinw.single_ion( ), summary: Stores single ion terms of the Hamiltonian.,
  keywords: sample, sidebar: sw_sidebar, permalink: spinw_single_ion.html, folder: spinw,
  mathjax: 'true'}

---
Stores single ion terms of the Hamiltonian.
Sub fields are:
  aniso   vector contains 1 x nMagAtom integers, each integer
          assignes one of the nMatrix from the .matrix field
          to a magnetic atom in the spinw.matom list as a single
          ion anisotropy (zeros for no anisotropy)
  g       vector contains 1 x nMagAtom integers, each integer
          assignes one of the nMatrix from the .matrix field
          to a magnetic atom in the spinw.matom list as a
          g-tensor
  field   external magnetic field stored in a 1x3 vector,
          default unit is Tesla
  T       temperature, scalar, default unit is Kelvin
 
See also SPINW.ADDANISO, SPINW.ADDG, SPINW.GETMATRIX, SPINW.SETMATRIX, SPINW.INTMATRIX.
