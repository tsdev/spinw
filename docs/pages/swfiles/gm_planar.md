---
{title: gm_planar, link: gm_planar, summary: planar magnetic structure constraint
    function, keywords: sample, sidebar: sw_sidebar, permalink: gm_planar, folder: swfiles,
  mathjax: 'true'}

---
  
### Syntax
  
`[s, k, n, name, pname, limit] = gm_planar(S0, x)`
  
### Description
  
`[s, k, n, name, pname, limit] = gm_planar(S0, x)` generates the
parameters of arbitrary planar magnetic structure from $$\varphi$$ angles
(in radian), ordering wave vector (rlu) and spin plane normal vector
($$xyz$$).
   
  
### Input Arguments
  
`x`
: Input parameters in the following order: 
  $$[\varphi_1, \varphi_2, ... , k_x, k_y, k_z, n_\theta, n_\varphi]$$.
  
`S0`
: Spin quantum number in a row vector $$(S_1, S_2, ...)$$ or scalar if all
  spins are equal.
  
### Output Arguments
  
`S`
: Matrix, containing the spin orientations with dimensions of $$[3\times n_{magExt}]$$.
      Every column contains the $$(S_x S_y S_z)$$ spin components of
      a magnetic atom in the $$xyz$$ coordinate system.
 
`k`
: Magnetic ordering wavevector in rlu units in a row vector.
 
`n`
: Normal vector around which the spins are rotating for non-zero
      propagation vector in a row vector.
 
`name`
: String, storing the name of the function.
 
`pname`
: Name of the input parameters in a cell: `{'Phi1_rad', ...}`.
 
`limit`
: Limits on the input parameters in a matrix with dimensions of $$[2\times n_{param}]$$. Every
      column contains a lower and upper limit on the corresponding
      parameter.
  
### See Also
  
[gm_spherical3d](gm_spherical3d) \| [gm_planard](gm_planard)
 
*[rlu]: Reciprocal Lattice Unit
 

{% include links.html %}
