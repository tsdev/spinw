---
{title: sw_rotmat( ), summary: rotates vectors around arbitrary axis in 3D, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_rotmat.html, folder: swfiles, mathjax: 'true'}

---
 
rotM = SW_ROTMAT(rotAxis, rotAngle)
 
It rotates vectors in V around rotAxis by rotAngle radian (positive angle
is the right-hand direction).
 
Input:
 
rotAxis   Axis of rotation, dimensions are [1 3].
rotAngle  Angle of rotation in radian (can be vector with dimensions of
          [1 nAng]).
 
Output:
 
rotM      Rotation matrix, dimensions are [3 3]. If rotAngle is a vector,
          rotM contains rotation matrices for every angle, it's
          dimensions are [3 3 nAng].
 
The rotation matrix defines rotations in a right-handed coordinate
system, the positive direction is counter-clockwise, when looking from
where the rotation axis points. To rotate any column vector use the
following:
  vp = rotM * v;
 
To rotate tensors (3x3 matrices) use the following command:
  Ap = rotM * A * rotM';
 
See also SPINW.GENMAGSTR, SW_ROT, SW_MIRROR.
 

