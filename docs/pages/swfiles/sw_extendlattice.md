---
{title: sw_extendlattice, link: sw_extendlattice, summary: creates superlattice, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_extendlattice, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`aList = sw_extendlattice(nExt,aList)`
  
`[aList, SSext] = sw_extendlattice(nExt,aList,SS)`
 
### Description
  
`aList = sw_extendlattice(nExt,aList)` creates a superlattice
and calculates all atomic positions within the new superlattice by
tiling it with the original cell.
 
`[aList, SSext] = sw_extendlattice(nExt,aList,SS)` also calculates the
bond matrix for the supercell by properly including all internal bonds
and bonds between atoms in different supercells.
  
### Input Arguments
  
`nExt`
: Size of the supercell in units of the original cell in a row vector
  with 3 elements.
  
`aList`
: List of the atoms, produced by [spinw.matom](spinw_matom).
  
`SS`
: Interactions matrices in the unit cell. Struct where each field
  contains an interaction matrix.
  
### Output Arguments
  
`aList`
: Parameters of the magnetic atoms in a struct with the following fields:
  * `RRext` Positions of magnetic atoms in lattice units of the supercell stored in a matrix with dimensions of $$[3\times n_{magExt}]$$.
  * `Sext`  Spin length of the magnetic atoms in a row vector with $$n_{magExt}$$ number of elements.
 
`SSext`
: Interaction matrix in the extended unit cell, struct type.
  In the struct every field is a matrix. Every column of the
  matrices describes a single bond, the following fields are generally
  defined:
	* `iso`     Isotropic exchange interactions.
	* `ani`     Anisotropic exchange interations.
	* `dm`      Dzyaloshinsky-Moriya interaction terms.
	* `gen`     General $$[3\times 3]$$ matrix contains the exchange interaction.
  
### See Also
  
[spinw.intmatrix](spinw_intmatrix)
 

{% include links.html %}
