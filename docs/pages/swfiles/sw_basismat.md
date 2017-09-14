---
{title: sw_basismat( ), link: sw_basismat, summary: determines allowed tensor components
    in a given point group symmetry, keywords: sample, sidebar: sw_sidebar, permalink: sw_basismat.html,
  folder: swfiles, mathjax: 'true'}

---
 
[M, asym] = SW_BASISMAT(symOp, r, tol) 
 
It determines the allowed matrix elements compatible with a given point
group symmetry. The matrix can describe exchange interaction or single
ion anisotropy.
 
Input:
symOp     Generators of the point group symmetry, dimensions are
          [3 3 nSym]. Each symOp(:,:,ii) matrix defines a rotation.
r         Distance vector between the two interacting atoms. For
          anisotropy r=0, dimensions are [3 1].
tol       Tolerance, optional. Default value is 1e-5.
 
Output:
 
M         Matrices, that span out the vector space of the symmetry
          allowed matrices, dimensions are [3 3 nM]. Any matrix is
          allowed that can be expressed as a linear combination of the
          symmetry allowed matrices.
asym      Logical vector, for each 3x3 matrix in M, tells whether it is
          antisymmetric, dimensions are [1 nM].
 
See also SPINW.GETMATRIX, SPINW.SETMATRIX.
 

