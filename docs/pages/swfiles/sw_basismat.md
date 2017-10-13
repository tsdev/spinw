---
{title: sw_basismat, link: sw_basismat, summary: determines allowed tensor components
    in a given point group symmetry, keywords: sample, sidebar: sw_sidebar, permalink: sw_basismat,
  folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`[m, asym] = sw_basismat(symop, r, tol)`
  
### Description
  
`[m, asym] = sw_basismat(symop, r, tol)` determines the allowed tensor
elements compatible with a given point group symmetry. The tensor can
describe exchange interaction or single ion anisotropy. The function
applies the symmetry invariance of the classical energy
$$\mathbf{S}_i\cdot \mathcal{M}\cdot \mathbf{S}_j$$. Thus this symmetry
analysis includes the transformation properties of spin operators as
well.
  
### Input Arguments
  
`symOp`
: Generators of the point group symmetry, in a matrix with dimensions of
  $$[3\times 3\times n_{sym}]$$ where each `symOp(:,:,ii)` matrix defines a rotation.
  
`r`
: Distance column vector between the two interacting atoms. For
  anisotropy $$r=0$$.
  
`tol`
: Tolerance, optional, default value is $$10^{-5}$$.
  
### Output Arguments
  
`M`
: Matrices, that span out the vector space of the symmetry
          allowed matrices, dimensions are $$[3\times 3\times n_M]$$. Any matrix is
          allowed that can be expressed as a linear combination of the
          symmetry allowed matrices.
 
`asym`
: Logical vector, for each $$[3\times 3]$$ matrix in $$M$$, tells whether it is
          antisymmetric stored in a row vector with $$n_M$$ elements.
  
### See Also
  
[spinw.getmatrix](spinw_getmatrix) \| [spinw.setmatrix](spinw_setmatrix)
 

{% include links.html %}
