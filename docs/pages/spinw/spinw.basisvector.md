---
{title: spinw.basisvector method, link: spinw.basisvector, summary: generates basis
    vectors of the crystal lattice, keywords: sample, sidebar: sw_sidebar, permalink: spinw_basisvector.html,
  folder: spinw, mathjax: 'true'}

---

### Syntax

`basisvec = basisvector(obj, {norm})`

### Description



### Examples

To change coordinate system:
relative atomic positions --> xyz
  r_xyz = basisvector * [ra; rb; rc];
reciprocal lattice units --> Å$$^{-1}$$ (xyz coordinate system)
  Q_xyz =  [h k l] * 2*pi*inv(basisvector);

### Input Arguments

`obj`
:pinw] object.

`norm`
:    If true, the basis vectors are normalized to 1, otherwise the
     length is equal to the lattice constants. Default is false.
     Optional.

### Output Arguments

basisVec  Stores the three basis vectors in columns, dimensions are 
          [3 3].
The 3x3 basisVec matrix can be used also as a coordinate transformation
matrix from the relative atomic position to positions in the xyz
coordinate system in Å units.

### See Also

[spinw](spinw.html) \| [spinw.abc](spinw_abc.html) \| [spinw.rl](spinw_rl.html)

{% include links.html %}
