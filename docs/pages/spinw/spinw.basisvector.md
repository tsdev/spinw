---
{title: spinw.basisvector method, link: spinw.basisvector, summary: generates lattice
    vectors, keywords: sample, sidebar: sw_sidebar, permalink: spinw_basisvector,
  folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`basisVec = basisvector(obj, {norm})`
  
### Description
  
`basisVec = basisvector(obj, {norm})` returns the lattice vectors of the
unit cell in a $$[3\times 3]$$ matrix, with the $$a$$, $$b$$ and $$c$$ vectors
stored in columns. The vectors are normalized to the lattice parameters
by default.
  
### Examples
  
The `basisVec` matrix can be used to change coordinate system, converting
between positions expressed in lattice units to positions expressed in
Å, using `r_lu` for lattice unit coordinates and `r_xyz` for
Å units (both stored in a column vector) the conversions are the
following:
```matlab
r_xyz = basisVec * r_lu
```
or
```matlab
r_lu = inv(basisVec)*r_xyz
```
 
It is also possible to convert between momentum vector in reciprocal
lattice units (rlu) into Å$$^{-1}$$ units. Assuming that momentum
vectors are row vectors:
```matlab
Q_xyz =  Q_rlu * 2*pi*inv(basisVec)
```
or
```matlab
Q_rlu = 1/(2*pi)*Q_xyz*basisVect
```
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
`norm`
: If `true`, the basis vectors are normalized to 1, otherwise the
  length of the basis vectors are equal to the lattice constants. Default
  is `false`.
  
### Output Arguments
  
`basisVec`
: Stores the three lattice vectors in columns, dimensions are $$[3\times 3]$$.
  
### See Also
  
[spinw](spinw) \| [spinw.abc](spinw_abc) \| [spinw.rl](spinw_rl)
 

{% include links.html %}
