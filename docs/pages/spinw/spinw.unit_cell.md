---
{title: spinw.unit_cell property, link: spinw.unit_cell, summary: Stores the atoms
    in the crystallographic unit cell., keywords: sample, sidebar: sw_sidebar, permalink: spinw_unit_cell.html,
  folder: spinw, mathjax: 'true'}

---
Sub fields are:
  r       positions of the atoms in the unit cell, in a
          3 x nAtom matrix, in lattice units
  S       spin quantum number of the atoms, in a 1 x nAtom
          vector, non-magnetic atoms have S=0
  label   label of the atom, strings in a 1 x nAtom cell
  color   color of the atom in 3 x nAtom matrix, where every
          column is an 0-255 RGB color
  ox      oxidation number of the atom, in a 1 x nAtom matrix
  occ     site occupancy in a 1 x nAtom matrix
  b       scattering length of the site for neutron and x-ray
          stored in a 2 x nAtom matrix, first row is neutron,
          second row is for x-ray
  ff      form factor of the site stored in a 2 x 9 x nAtom
          matrix, first row is the magnetic form factor for
          neutrons, the second row is the charge form factor
          for x-ray cross section
  Z       atomic number
  A       atomic mass (N+Z) for isotopes and -1 for natural
          distribution of isotopes
  biso    Isotropic displacement factors in units of Angstrom^2.
          Definition is the same as in FullProf, defining the
          Debye-Waller factor as:
              Wd = 1/8*biso/d^2
          including in the structure factor as exp(-2Wd)
 
See also SPINW.ADDATOM, SPINW.ATOM, SPINW.MATOM, SPINW.NEWCELL, SPINW.PLOT.

