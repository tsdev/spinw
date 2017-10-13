---
{title: sw_rot, link: sw_rot, summary: rotates vectorsin 3D, keywords: sample, sidebar: sw_sidebar,
  permalink: sw_rot, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`[~, R] = sw_rot(rotAxis,rotAngle)`
  
`[VR, R] = sw_rot(rotAxis,rotAngle,V)`
 
### Description
  
`[~, R] = sw_rot(rotAxis,rotAngle)` produces the `R` rotation matrix that
rotates any vector around the given `rotAxis` rotation axis by `rotAngle`
angle in radian. Positive rotation is the right-hand direction around the
rotation axis and using the following rotation formula:
```matlab
VR = R*V
```
 
To rotate tensors ($$3\times 3$$ matrices) use the following formula:
```matlab
Mp = R * M * R';
```
 
`[VR, R] = sw_rot(rotAxis,rotAngle,V)` also rotates the given `V`
vectors where `VR` are the transformed vectors.
  
### Input Arguments
  
`rotAxis`
: Axis of rotation, stored in a row vector with 3 elements.
  
`rotAngle`
: Angle of rotation in radian, can be also a row vector with $$n_{ang}$$
  number of elements.
  
`V`
: Matrix with 3 rows, where each column is a vector in 3D space.
  
### Output Arguments
  
`VR`
: Rotated vectors, stored in a matrix with dimensions of $$[3\times N
  n_{ang}]$$.
 
`R`
: Rotation matrix with dimensions of $$[3\times 3]$$ if a single rotation
  angle is given. If `rotAngle` is a vector, `R` will contain a
  rotation matrix for each angle, it's dimensions are $$[3\times 3\times
  n_{ang}]$$.
 
### See Also
  
[sw_rotmat](sw_rotmat) \| [sw_mirror](sw_mirror)
 

{% include links.html %}
