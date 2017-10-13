---
{title: swsym.bond, link: swsym.bond, summary: generates all symmetry equivalent bonds,
  keywords: sample, sidebar: sw_sidebar, permalink: swsym_bond, folder: swsym, mathjax: 'true'}

---
  
### Syntax
  
`[genBond, uBond] = swsym.bond(r,bv,bond,symOp)`
  
`[genBond, uBond] = swsym.bond(r,bv,bond,symOp,tol)`
 
### Description
  
`[genBond, uBond] = swsym.bond(r,bv,bond,symOp)` generates all bonds that
are symmetry equivalent to the given `bond`. The function uses the given
space group operators and positions of magnetic atoms to return a list of
equivalent bonds in a matrix. The function also checks the validity of
the calculation by measuring the length of each equivalent bond using the
given `bv` base and if the difference in length between equivalent bonds
is larger than the tolerance throws a warning.
  
`[genBond, uBond] = swsym.bond(r,bv,bond,symOp,tol)` also defines the
tolerance using `tol`.
 
### Input Arguments
  
`r`
: Positions of the magnetic atoms in lattice units stored in a matrix
  with dimensions of $$[3\times n_{magAtom}]$$.
  
`bv`
: Basis vectors that define the lattice, used for checking the bond
  length of equivalent bonds, see [spinw.basisvector](spinw_basisvector) for details.
  
`bond`
: Vector that contains the starting bond with elements of 
  `[dl_a dl_b dl_c atom_1 atom_2]`, where `dl` is vector of lattice
  translation between the two atoms if they are not in the same unit cell
  in lattice units, `atom_1` and `atom_2` are indices of atoms in the
  list of positions stored in `r`.
  
`symOp`
: Matrix, that contains the rotation and translation operators of
  the space group with dimensions of $$[3\times 4\times n_{op}]$$.
  
`tol`
: Tolerance, default value is $$10^{-5}$$.
  
### Output Arguments
  
`genBond`
: Matrix, whith each column defines a bond, the meaning of each
          row is the same as the input `bond` variable.
 
`uBond`
: Logical variable, `true` if all the generated bonds are unique.
  
### See Also
  
[spinw.gencoupling](spinw_gencoupling) \| [swsym.operator](swsym_operator) \| [swsym.position](swsym_position)
 

{% include links.html %}
