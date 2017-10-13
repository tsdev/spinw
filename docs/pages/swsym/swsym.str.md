---
{title: swsym.str, link: swsym.str, summary: generates a string equivalent of symmetry
    operators, keywords: sample, sidebar: sw_sidebar, permalink: swsym_str, folder: swsym,
  mathjax: 'true'}

---
  
### Syntax
  
`symstr = swsym.str(symop)`
  
### Description
  
`symStr = swsym.str(symOp)` generates a string equivalent of the given
symmetry operator matrix. The string contains the operators separated by
`;` and the $$xyz$$ axis transformations are separated by `,`. For example
a valid symmetry operator is `'x,y+1/2,z+1/2'`. Translations are given as
fractions and the `'xyz'` letters correspond to the 3 crystal axes.
  
### Input Arguments
  
`symOp`
: Symmetry operators in a matrix with dimensions of $$[3\times 4\times
  n_{op}]$$, where rotations matrices are stored in `symOp(:,1:3,:)` and
  and translation vectors in `symOp(:,4,:)`.
  
### Output Arguments
  
`strSym`
: String that contains the symmetry operations.
  
### See Also
  
[swsym.add](swsym_add) \| [swsym.generator](swsym_generator)
 

{% include links.html %}
