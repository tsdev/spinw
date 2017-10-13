---
{title: swsym package, link: swsym package, summary: package to handle symmetry operations,
  keywords: sample, sidebar: sw_sidebar, permalink: swsym, folder: swsym, mathjax: 'true'}

---
 
This package deals with symmetry operators of crystallographic space
groups. It can read the standard space group definitions stored in
`symmetry.dat`, generate all symmmetry elements, determine all symmetry
equivalent positions, etc. 
 
All symmetry operators `symOp` are defined by a matrix with dimensions of
$$[3\times 4\times n_{op}]$$, where `symOp(1:3,1:3,:)` stores the $$[3\times
3]$$ rotation matrices while the `symOp(1:3,4,:)` holds the corresponding
translation vectors.
 
### Files
 
* [swsym.add](swsym_add) saves user defined symmetry operators
* [swsym.bond](swsym_bond) generates all symmetry equivalent bonds
* [swsym.generator](swsym_generator) returns symmetry operators of a given space group
* [swsym.genreduce](swsym_genreduce) reduces symmetry operators to the generators
* [swsym.isop](swsym_isop) determines if a matrix is symmetry operator
* [swsym.operator](swsym_operator) generates all symmetry elements from given space group
* [swsym.oporder](swsym_oporder) determine the order of the symmetry operator
* [swsym.point](swsym_point) determines local point group symmetry in a space group
* [swsym.position](swsym_position) generates symmetry equivalent positions
* [swsym.str](swsym_str) generates a string equivalent of symmetry operators

{% include links.html %}
