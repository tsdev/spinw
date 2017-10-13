---
{title: spinw.quickham method, link: spinw.quickham, summary: quickly generate magnetic
    Hamiltonian, keywords: sample, sidebar: sw_sidebar, permalink: spinw_quickham,
  folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`quickham(obj,J)`
  
### Description
  
`quickham(obj,J)` generates the bonds from the predefined crystal
structure and assigns exchange values to bonds such as `J(1)` to first
neighbor, `J(2)` for second neighbor etc. The command will erase all
previous bonds, anisotropy, g-tensor and matrix definitions. Even if
`J(idx) == 0`, the corresponding bond and matrix will be created.
   
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
`J`
: Vector that contains the Heisenberg exchange values. `J(1)` for
     first neighbor bonds, etc.
 
### See Also
 
[spinw.gencoupling](spinw_gencoupling) \| [spinw.addcoupling](spinw_addcoupling) \| [spinw.matrix](spinw_matrix) \|
[spinw.addmatrix](spinw_addmatrix)
 

{% include links.html %}
