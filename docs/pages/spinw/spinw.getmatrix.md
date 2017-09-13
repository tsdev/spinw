---
{title: spinw.getmatrix( ), summary: gives the symmetry allowed matrices for a given
    coupling or anisotropy, keywords: sample, sidebar: sw_sidebar, permalink: spinw_getmatrix.html,
  folder: spinw, mathjax: 'true'}

---
gives the symmetry allowed matrices for a given coupling or anisotropy
 
aMat = GETMATRIX(obj, 'Option1', Value1, ...)
 
Input:
 
obj           spinw class object.
 
Options:
 
One of the following options has to be given in the input:
 
mat           Label or index of the matrix that is already assigned to
              a bond, anisotropy or g-tensor.
bond          Index of the bond in spinw.coupling.idx.
subIdx        Selects a certain bond, within efault value is 1.
aniso         Label or index of the magnetic atom that has a single ion
              anisotropy matrix is assigned.
gtensor       Label or index of the magnetic atom that has a g-tensor is 
              assigned.
 
Optional inputs:
 
tol       Tolerance for printing the output for the smallest matrix
          element.    
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
 
aMat      If no prefactors are defined, aMat contains all symmetry
          allowed elements of the coupling/anisotropy matrix, dimensions
          are [3 3 nSymMat]. If prefactor is defined, it is a single 3x3
          matrix, that is a sum of all symmetry allowed elemenets
          multiplied by the given prefactors.
 
Example:
 
cryst = spinw;
cryst.genlattice('sym','P 4')
cryst.addatom('r',[0 0 0],'label','MCu2')
cryst.addmatrix('label','A','value',eye(3))
cryst.gencoupling
cryst.addaniso('A')
cryst.getmatrix('mat','A');
 
The above example determines the allowed anisotropy matrix elements in
the C4 point group symmetry (the symmetry at the [0 0 0] atomic
position) and prints them onto the Command Window. The allowed matrix
elements are: diag([A A B]).
 
See also SPINW.SETMATRIX.
 
