---
{title: Package swsym, summary: The SWSYM library handles the symmetry operations
    for SpinW objects. The, keywords: sample, sidebar: sw_sidebar, permalink: swsym.html,
  folder: +swsym, mathjax: 'true'}

---
The SWSYM library handles the symmetry operations for SpinW objects. The
functions can be used independently from other parts of SpinW. All
symmetry operators (symOp) are defined by a matrix with dimensions [3 4
nOp], where symOp(:,1:3,:) defines the rotation matrices while the
symOp(:,4,:) the corresponding translations. Also the standard settings
of the space groups are stored in the symmetry.dat file that can be
loaded using the SWSYM.generator or SWSYM.operator functions.
 
### Files

* [swsym.add](/swsym_add) saves user defined symmetry operators
* [swsym.bond](/swsym_bond) generates all symmetry equivalent bonds
* [swsym.generator](/swsym_generator) returns symmetry operators of a given space group
* [swsym.genreduce](/swsym_genreduce) reduces the list of symmetry operators to the generators
* [swsym.isop](/swsym_isop) function determines whether the matrix is a symmetry operator
* [swsym.operator](/swsym_operator) calculates all symmetry operators or general positions for a space group
* [swsym.oporder](/swsym_oporder) determine the order of the symmetry operator
* [swsym.point](/swsym_point) determines point group symmetry at a given position
* [swsym.position](/swsym_position) generates symmetry equivalent positions
* [swsym.str](/swsym_str) generates a string equivalent of symmetry operators
