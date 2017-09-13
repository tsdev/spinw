---
{title: spinw.setmatrix( ), summary: changes the selected matrix inside the spinw
    object, keywords: sample, sidebar: sw_sidebar, permalink: spinw_setmatrix.html,
  folder: spinw, mathjax: 'true'}

---
changes the selected matrix inside the spinw object
 
SETMATRIX(obj, 'Option1', Value1, ...)
 
Input:
 
obj           spinw class object.
 
Options:
 
One of the below options has to be given:
 
mat           Label or index of the matrix that is already assigned to
              a bond, anisotropy or g-tensor.
bond          Index of the bond in spinw.coupling.idx.
subIdx        Selects a certain bond, within efault value is 1.
aniso         Label or index of the magnetic atom that has a single ion
              anisotropy matrix is assigned.
gtensor       Label or index of the magnetic atom that has a g-tensor is 
              assigned.
 
Optional inputs:
 
pref      Defines prefactors as a vector for the symmetry allowed
          components, dimensions are [1 nSymMat]. Alternatively, if only
          a few of the symmetry allowed matrices have non-zero
          prefactors, use:
              {[6 0.1 5 0.25]}
          This means, the 6th symmetry allowed matrix have prefactor 0.1,
          the 5th symmetry allowed matrix have prefactor 0.25. Since
          Heisenberg isotropic couplings are always allowed, a cell with
          a single element will create a Heisenberg coupling, example:
              {0.1}
          This is identical to obj.matrix.mat = eye(3)*0.1
          For DM interactions (antisymmetric coupling matrices), use
          three element vector in the cell:
              {[D1 D2 D3]}
          In this case, these will be the prefactors of the 3
          antisymmetric symmetry allowed matrices. In case no crystal
          symmetry is defined, these will define directly the components
          of the  DM interaction in the xyz coordinate system. Be
          carefull with the sign of the DM interaction, it depends on the
          order of the two interacting atoms! Default value is {1}.
          For anisotropy matrices antisymmetric matrices are not allowed.
 
Output:
 
The selected obj.matrix.mat will contain the new value.
 
Example:
 
...
setmatrix(crystal,'label','J1','pref',{[6 0.235]})
 
This example will set 'J1' coupling to the 6th symmetry allowed matrix,
with prefactor 0.235.
 
setmatrix(crystal,'label','J2','pref',{1.25})
 
This will set 'J2' to antiferromagnetic Heisenberg exchange, with value
of 1.25 meV.
 
See also SPINW, SPINW.GENCOUPLING, SPINW.GETMATRIX.
 
