---
{title: spinw.rl( ), summary: generates reciprocal lattice basis vectors of the crystal
    lattice, keywords: sample, sidebar: sw_sidebar, permalink: spinw_rl.html, folder: spinw,
  mathjax: 'true'}

---
 
rlVec = RL(obj, {norm})
 
Input:
 
obj       spinw class object.
norm      If true, the basis vectors are normalized to 1. Default is false.
          Optional.
 
Output:
 
rlVec     Stores the three basis vectors in columns, dimensions are
          [3 3].
 
The 3x3 rlVec matrix can be used also as a coordinate transformation
matrix from the relative atomic position to positions in the xyz
coordinate system in Angstrom units.
 
Example:
 
To convert from reciprocal lattice unit to Angstrom^-1 (xyz coordinate system):
 
  Q_xyz =  [h k l] * rlVect;
 
See also SPINW, SPINW.ABC, SPINW.BASISVECTOR.
 

