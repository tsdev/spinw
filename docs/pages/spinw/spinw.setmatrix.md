---
{title: spinw.setmatrix method, link: spinw.setmatrix, summary: changes the selected
    matrix inside the spinw object, keywords: sample, sidebar: sw_sidebar, permalink: spinw_setmatrix.html,
  folder: spinw, mathjax: 'true'}

---

### Syntax

`setmatrix(obj,Name,Value)`

### Description



### Examples

...
setmatrix(crystal,'label','J1','pref',{[6 0.235]})
This example will set 'J1' coupling to the 6th symmetry allowed matrix,
with prefactor 0.235.
setmatrix(crystal,'label','J2','pref',{1.25})
This will set 'J2' to antiferromagnetic Heisenberg exchange, with value
of 1.25 meV.

### Input Arguments

`obj`
: [spinw](spinw.html) object.

### Name-Value Pair Arguments

`'One'`
:below options has to be given:

`'mat'`
:   Label or index of the matrix that is already assigned to
    a bond, anisotropy or g-tensor.

`'bond'`
:   Index of the bond in spinw.coupling.idx.

`'subIdx'`
:   Selects a certain bond, within efault value is 1.

`'aniso'`
:   Label or index of the magnetic atom that has a single ion
    anisotropy matrix is assigned.

`'gtensor'`
:   Label or index of the magnetic atom that has a g-tensor is 
    assigned.

`'Optional'`
:puts:

`'pref'`
:efines prefactors as a vector for the symmetry allowed
 omponents, dimensions are [1 nSymMat]. Alternatively, if only
  few of the symmetry allowed matrices have non-zero
 refactors, use:
    {[6 0.1 5 0.25]}
 his means, the 6th symmetry allowed matrix have prefactor 0.1,
 he 5th symmetry allowed matrix have prefactor 0.25. Since
 eisenberg isotropic couplings are always allowed, a cell with
  single element will create a Heisenberg coupling, example:
    {0.1}
 his is identical to obj.matrix.mat = eye(3)*0.1
 or DM interactions (antisymmetric coupling matrices), use
 hree element vector in the cell:
    {[D1 D2 D3]}
 n this case, these will be the prefactors of the 3
 ntisymmetric symmetry allowed matrices. In case no crystal
 ymmetry is defined, these will define directly the components
 f the  DM interaction in the xyz coordinate system. Be
 arefull with the sign of the DM interaction, it depends on the
 rder of the two interacting atoms! Default value is {1}.
 or anisotropy matrices antisymmetric matrices are not allowed.

### Output Arguments

The selected obj.matrix.mat will contain the new value.

### See Also

[spinw](spinw.html) \| [spinw.gencoupling](spinw_gencoupling.html) \| [spinw.getmatrix](spinw_getmatrix.html)

