---
{title: swsym.isop, link: swsym.isop, summary: determines if a matrix is symmetry
    operator, keywords: sample, sidebar: sw_sidebar, permalink: swsym_isop, folder: swsym,
  mathjax: 'true'}

---
  
### Syntax
  
`result = swsym.isop(op)`
  
### Description
  
`result = swsym.isop(op)` determines whether the given matrix has
dimensions that is compatible with the size requirements of space group
operators. The given `op` matrix has to have dimensions of $$[3\times
4\times n_{op}]$$. The function returns `true` only if the input has these
dimensions.
 
### See Also
  
[swsym.generator](swsym_generator) \| [swsym.operator](swsym_operator)
 

{% include links.html %}
