---
{title: gm_spherical3d, link: gm_spherical3d, summary: magnetic structure constraint
    function with spherical parameterisation, keywords: sample, sidebar: sw_sidebar,
  permalink: gm_spherical3d, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`[m, k, n, name, pname, limit] = gm_spherical3d(S0, x)`
  
### Description
  
`[m, k, n, name, pname, limit] = gm_spherical3d(S0, x)` generates
magnetic structure from given parameters while constraining the length of
the spin on each atom. The parametrization of the magnetic structure
consists of 2 spherical coordinates $$(\theta,\varphi)$$ angles per
magnetic atom. All angles are in radian.
    
### Input Arguments
    
`x`
: Input parameters in the following order: 
  $$[\theta_1, \varphi_1, ... , k_x, k_y, k_z, n_\theta, n_\varphi]$$.
  
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
    
[gm_planar](gm_planar)

{% include links.html %}
