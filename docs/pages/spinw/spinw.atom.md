---
{title: spinw.atom method, link: spinw.atom, summary: generates all atomic positions
    in the unit cell, keywords: sample, sidebar: sw_sidebar, permalink: spinw_atom.html,
  folder: spinw, mathjax: 'true'}

---

### Syntax

`atomlist = atom(obj)`

### Description

From the given atomic positions in the obj.unit_cell field the function
generates all atomic positions using the symmetry operators of the given
space group. If no symmetry is defined, the function returns simply the
positions stored in obj.unit_cell.
 

### Examples

cryst = spinw;
cryst.genlattice('lat_const',[8 8 8],'spgr','x+1/2,y+1/2,z;x+1/2,y,z+1/2;x,y+1/2,z+1/2','label','FCC')
cryst.addatom('r',[0 0 0],'label','Atom1')
atomList = cryst.atom;
This will create a new space group, that contains all the translations of
the FCC lattice. Then creates a crystal with an atom at [0 0 0] position.
The cryst.atom lists all 4 symmetry equivalent positions generated using
the 'FCC' symmetry operators.

### Input Arguments

`obj`
: [spinw](spinw.html) object.

### Output Arguments

atomList is a structure with the following fields:
  r       Positions of the atoms in lattice units, dimensions are 
          [3 nAtom]. 
  idx     Pointer to the atom in the unit_cell field, dimensions are
          [nAtom 1].
  mag     Logical variable, whether the spin of the atom is non-zero,
          dimensions are [nAtom 1].

### See Also

[spinw](spinw.html) \| [spinw.matom](spinw_matom.html) \| [swsym.add](swsym_add.html) \| [spinw.genlattice](spinw_genlattice.html) \| [spinw.addatom](spinw_addatom.html)

