---
{title: spinw.mag_str( ), summary: Stores the magnetic structure., keywords: sample,
  sidebar: sw_sidebar, permalink: spinw_mag_str.html, folder: spinw, mathjax: 'true'}

---
Stores the magnetic structure.
Sub fields are:
  S       stores the moment direction for every spin in the
          crystallographic or magnetic supercell in a
          3 x nMagExt matrix, where nMagExt = nMagAtom*prod(N_ext)
  k       magnetic ordering wave vector in a 3x1 vector
  n       normal vector to the rotation of the moments in
          case of non-zero ordering wave vector, dimensions are
          3x1
  N_ext   Size of the magnetic supercell, default is [1 1 1]
          if the magnetic cell is identical to the
          crystallographic cell, the 1x3 vector extends the
          cell along the a, b and c axis
 
See also SPINW.GENMAGSTR, SPINW.OPTMAGSTR, SPINW.ANNEAL, SPINW.MOMENT, SPINW.NMAGEXT, SPINW.STRUCTFACT.

