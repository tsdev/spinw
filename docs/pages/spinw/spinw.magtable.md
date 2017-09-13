---
{title: spinw.magtable( ), summary: creates tabulated list of all magnetic moments
    stored in obj, keywords: sample, sidebar: sw_sidebar, permalink: spinw_magtable.html,
  folder: spinw, mathjax: 'true'}

---
creates tabulated list of all magnetic moments stored in obj
 
moments = MAGTABLE(obj)
 
The function lists the APPROXIMATED moment directions (using the rotating
coordinate system notation) in the magnetic supercell, whose size is
defined by the obj.mag_str.nExt field. The positions of the magnetic
atoms are in lattice units.
 
Input:
 
obj           spinw class object.
 
Output:
 
'moments' is struct type data that contains the following fields:
  M           Matrix, where every column defines a magnetic moment,
              dimensions are [3 nMagExt].
  e1,e2,e3    Unit vectors of the coordinate system used for the spin
              wave calculation, the i-th column contains a normalized
              vector for the i-th moment. e3 is parallel to the magnetic
              moment, e1 and e2 span a right handed orthogonal coordinate
              system.
  R           Matrix, where every column defines the position of the
              magnetic atom in lattice units.
  atom        Pointer to the magnetic atom in the subfields of
              spinw.unit_cell.
 
See also SPINW.GENMAGSTR.
 
