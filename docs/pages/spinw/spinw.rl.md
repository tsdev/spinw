---
{title: spinw.rl method, link: spinw.rl, summary: generates reciprocal lattice vectors,
  keywords: sample, sidebar: sw_sidebar, permalink: spinw_rl, folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`rlVec = rl(obj, {norm})`
  
### Description
  
`rlVec = rl(obj, {norm})` returns the lattice vectors of the reciprocal
lattice in a $$[3\times 3]$$ matrix, with the $$a^*$$, $$b^*$$ and $$c^*$$ vectors
stored in **rows**. 
 
  
### Examples
  
To convert from reciprocal lattice unit to Å$$^{-1}$$ ($$xyz$$
Cartesian coordinate system) use: (`Q_rlu` is a row vector with 3
elements):
 
```matlab
Q_xyz =  Q_rlu * rlVect
```
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
`norm`
: If `true`, the basis vectors are normalized to 1. Default values is
`false`, optional.
  
### Output Arguments
  
`rlVec`
: Stores the three basis vectors in the rows of matrix with dimensions of
  $$[3\times 3]$$.
  
### See Also
  
[spinw](spinw) \| [spinw.abc](spinw_abc) \| [spinw.basisvector](spinw_basisvector)
 

{% include links.html %}
