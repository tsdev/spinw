---
{title: '@spinw.lattice( )', summary: Stores the unit cell parameters., keywords: sample,
  sidebar: sw_sidebar, permalink: spinw_lattice.html, folder: '@spinw', mathjax: 'true'}

---
Stores the unit cell parameters.
Sub fields are:
  lat_const   lattice constants in a 1x3 vector in Angstrom units
  angle       (alpha,beta,gamma) angles in 1x3 vector in radian
  sym         crystal space group, line number in symmetry.dat file
  origin      Origin of the cell in l.u.
  label       Label of the space group.
 
See also SPINW.GENLATTICE, SPINW.ABC, SPINW.BASISVECTOR, SPINW.NOSYM.
