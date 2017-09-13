---
{title: swsym.oporder( ), summary: determine the order of the symmetry operator, keywords: sample,
  sidebar: sw_sidebar, permalink: swsym_oporder.html, folder: swsym, mathjax: 'true'}

---
determine the order of the symmetry operator
 
N = SWSYM.OPORDER(symOp)
 
It determines the order of the symOp symmetry operator, where
symOp(:,1:3) is a rotation matrix and symOp(:,4) is a translation.
Maximum order is 10 if the matrix is not a rotation matrix of any
crystallographic point group.
 
Input:
 
symOp 	Symmetry operator in a matrix.
 
Example:
 
R^sw_symorder([R zeros(3,1)]) == eye(3);
 
See also SWSYM.GENERATOR, SW_BASISMAT.
 
