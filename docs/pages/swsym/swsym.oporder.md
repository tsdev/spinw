---
{title: swsym.oporder, link: swsym.oporder, summary: determine the order of the symmetry
    operator, keywords: sample, sidebar: sw_sidebar, permalink: swsym_oporder.html,
  folder: swsym, mathjax: 'true'}

---

### Syntax

` `

### Description

symOp(:,1:3) is a rotation matrix and symOp(:,4) is a translation.
Maximum order is 10 if the matrix is not a rotation matrix of any
crystallographic point group.
 

### Examples

R^sw_symorder([R zeros(3,1)]) == eye(3);

### Input Arguments

% `symOp`
:	Symmetry operator in a matrix.

### See Also

[swsym.generator](swsym_generator.html) and [sw_basismat](sw_basismat.html)

