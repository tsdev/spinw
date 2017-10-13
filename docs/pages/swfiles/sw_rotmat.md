---
{title: sw_rotmat, link: sw_rotmat, summary: generates 3D rotation matrix, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_rotmat, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`R = sw_rotmat(rotAxis,rotAngle)`
 
### Description
  
`R = sw_rotmat(rotAxis,rotAngle)` produces the `R` rotation matrix that
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
 
### Input Arguments
  
`rotAxis`
: Axis of rotation, stored in a row vector with 3 elements.
  
`rotAngle`
: Angle of rotation in radian, can be also a row vector with $$n_{ang}$$
  number of elements.
  
### Output Arguments
  
`R`
: Rotation matrix with dimensions of $$[3\times 3]$$ if a single rotation
  angle is given. If `rotAngle` is a vector, `R` will contain a
  rotation matrix for each angle, it's dimensions are $$[3\times 3\times
  n_{ang}]$$.
 
### See Also
  
[sw_rot](sw_rot) \| [sw_mirror](sw_mirror)
 

{% include links.html %}
