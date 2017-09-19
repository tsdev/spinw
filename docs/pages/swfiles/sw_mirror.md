---
{title: sw_mirror( ), link: sw_mirror, summary: mirrors a 3D vector, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_mirror.html, folder: swfiles, mathjax: 'true'}

---

### Syntax

` `

### Description

normal vector mNorm.
 

### Input Arguments

% `n`
: Vector, normal to the mirror plane, dimensions are [1 3].

% `V`
: Matrix of 3D vectors, dimensions are [3 N], optional.

### Output Arguments

V         Mirrored vectors, dimensions are [3 N].
mirM      Matrix of the mirror operation, dimensions are [3 3].
To mirror any column vector use the following:
vp = mirM * v;
To apply mirror plane operation on tensors (3x3 matrices) use the
following command:
Ap = mirM * A * mirM';

### See Also

[spinw.genmagstr](spinw_genmagstr.html) and [sw_rot](sw_rot.html)

