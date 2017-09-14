---
{title: swsym package, link: swsym package, summary: package to handle symmetry operations,
  keywords: sample, sidebar: sw_sidebar, permalink: swsym.html, folder: swsym, mathjax: 'true'}

---
 
The functions can be used independently from other parts of spinw. All
symmetry operators `symOp` are defined by a matrix with dimensions of
$$[3\times 4\times n_{op}]$$, where `symOp(:,1:3,:)` defines the rotation
matrices while the `symOp(:,4,:)` the corresponding translations. Also
the standard settings of the space groups are stored in the
`symmetry.dat` file that can be loaded using the [swsym.generator](swsym_generator.html) or
[swsym.operator](swsym_operator.html) functions.
 
### Files
 
* [swsym.add](swsym_add.html) saves user defined symmetry operators
* [swsym.bond](swsym_bond.html) generates all symmetry equivalent bonds
* [swsym.generator](swsym_generator.html) returns symmetry operators of a given space group
* [swsym.genreduce](swsym_genreduce.html) reduces the list of symmetry operators to the generators
* [swsym.isop](swsym_isop.html) function determines whether the matrix is a symmetry operator
* [swsym.operator](swsym_operator.html) calculates all symmetry operators or general positions for a space group
* [swsym.oporder](swsym_oporder.html) determine the order of the symmetry operator
* [swsym.point](swsym_point.html) determines point group symmetry at a given position
* [swsym.position](swsym_position.html) generates symmetry equivalent positions
* [swsym.str](swsym_str.html) generates a string equivalent of symmetry operators

