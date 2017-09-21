---
{title: spinw.rl method, link: spinw.rl, summary: generates reciprocal lattice basis
    vectors of the crystal lattice, keywords: sample, sidebar: sw_sidebar, permalink: spinw_rl.html,
  folder: spinw, mathjax: 'true'}

---

### Syntax

`rlvec = rl(obj, {norm})`

### Description



### Examples

To convert from reciprocal lattice unit to Å$$^{-1}$$ (xyz coordinate system):
  Q_xyz =  [h k l] * rlVect;

### Input Arguments

`obj`
:pinw] object.

`norm`
:    If true, the basis vectors are normalized to 1. Default is false.
     Optional.

### Output Arguments

rlVec     Stores the three basis vectors in columns, dimensions are
          [3 3].
The 3x3 rlVec matrix can be used also as a coordinate transformation
matrix from the relative atomic position to positions in the xyz
coordinate system in Å units.

### See Also

[spinw](spinw.html) \| [spinw.abc](spinw_abc.html) \| [spinw.basisvector](spinw_basisvector.html)

{% include links.html %}
