---
{title: swsym.bond, link: swsym.bond, summary: generates all symmetry equivalent bonds,
  keywords: sample, sidebar: sw_sidebar, permalink: swsym_bond.html, folder: swsym,
  mathjax: 'true'}

---

### Syntax

`[gencp ugencp] = swsym.bond(r, bv, bond, symop, {tol})`

### Description



### Input Arguments

`r`
: Positions of the magnetic atoms in the unit cell in a matrix
  with dimensions of [3 nMagAtom], in lattice units.

`bv`
: Basis vectors of the lattice, used for error checking of the
  calculation.

`bond`
: Vector, contains the starting bond with elements of 
  [dl_x dl_y dl_z atom_1 atom_2], where dl is vector of lattice
  translation between the two atom if they are not in the same
  cell, atom_1 and atom_2 are indices of atoms in the list of
  positions stored in r.

`symOp`
: Matrix, that contains the rotation and translation operators of
  the space group with dimensions of [3 4 nOp].

`tol`
: Tolerance, optional. Default value is 1e-5.

### Output Arguments

genCp     Matrix, each column contains a coupling, the meaning of each
          row are the same as the input coupling variable.
ugenCp    Logical variable, true if the coupling is unique in the list of
          generated couplings.

### See Also

[spinw.gencoupling](spinw_gencoupling.html) \| [swsym.operator](swsym_operator.html) \| [swsym.position](swsym_position.html)

