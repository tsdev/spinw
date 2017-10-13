---
{title: spinw.getmatrix method, link: spinw.getmatrix, summary: determines the symmetry
    allowed tensor elements, keywords: sample, sidebar: sw_sidebar, permalink: spinw_getmatrix,
  folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`amat = getmatrix(obj,Name,Value)`
  
### Description
  
`amat = getmatrix(obj,Name,Value)` determines the symmetry allowed
elements of the exchange, single-ion anistropy and g-tensor. For bonds,
the code first determines the point group symmetry on the center of the
bond and calculates the allowed eelements of the exchange tensor
accordingly. For anisotropy and g-tensor, the point group symmetry of the
selected atom is considered. For example the code can correctly calculate
the allowed Dzyaloshinskii-Moriya vectors.
  
### Examples
  
To following code will determine the allowed anisotropy matrix elements
in the $$C4$$ point group (the symmetry at the $$(0,0,0)$$ atomic position).
The allowed matrix elements will be `diag([A A B])`:
 
```matlab
cryst = spinw
cryst.genlattice('sym','P 4')
cryst.addatom('r',[0 0 0],'label','MCu2')
cryst.addmatrix('label','A','value',1)
cryst.gencoupling
cryst.addaniso('A')
cryst.getmatrix('mat','A')
```
*Output*
```
The symmetry analysis of the anisotropy matrix of atom 1 ('MCu2'):
 position (in lattice units): [0.000,0.000,0.000]
 label of the assigned matrix: 'A'
 allowed elements in the symmetric matrix:
  S = | A| 0| 0|
      | 0| A| 0|
      | 0| 0| B|
```
 
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
### Name-Value Pair Arguments
 
At least one of the following option has to be defined:
  
`mat`
: Label or index of a matrix that is already assigned to
  a bond, anisotropy or g-tensor, e.g. `J1`.
  
`bond`
: Index of the bond in `spinw.coupling.idx`, e.g. 1 for first neighbor
  bonds.
  
`aniso`
: Label or index of the magnetic atom that has a single ion
  anisotropy matrix is assigned, e.g. `Cr1` if `Cr1` is a magnetic atom.
  
`gtensor`
: Label or index of the magnetic atom that has a g-tensor is 
  assigned.
 
Optional inputs:
  
`subIdx`
: Selects a certain bond, within equivalent bonds. Default value is 1.
 
`tol`
: Tolerance for printing the output for the smallest matrix
  element.
  
`pref`
: If defined `amat` will contain a single $$[3\times 3]$$ matrix by
  multuplying the calculated tensor components with the given prefactors.
  Thus `pref` should contain the same number of elements as the number of
  symmetry allowed tensor components. Alternatively, if only a few of the
  symmetry allowed matrices have non-zero prefactors, use e.g. 
  `{[6 0.1 5 0.25]}` which means, the 6th symmetry allowed matrix have
  prefactor 0.1, the 5th symmetry allowed matrix have prefactor 0.25.
  Since Heisenberg isotropic couplings are always allowed, a cell with a
  single element will create a Heisenberg coupling, e.g. `{0.1}`, which is
  identical to `obj.matrix.mat = eye(3)*0.1`. For Dzyaloshinskii-Moriya
  interactions (antisymmetric exchange matrices), use a three element
  vector in a cell, e.g. `pref = {[D1 D2 D3]}`. In this case, these will
  be the prefactors of the 3 antisymmetric allowed matrices. In
  case no crystal symmetry is defined, these will define directly the
  components of the  Dzyaloshinskii-Moriya interaction in the xyz
  coordinate system.
 
  {% include note.html content=" Be carefull with the sign of the Dzyaloshinskii-Moriya
  interaction, it depends on the counting order of the two interacting
  atoms! For single-ion anisotropy and g-tensor antisymmetric matrices
  are forbidden in any symmetry." %}
  
### Output Arguments
  
`aMat`
: If no prefactors are defined, `aMat` contains all symmetry
  allowed elements of the selected tensor, dimensions are $$[3\times 3\times n_{symmat}]$$.
  If a prefactor is defined, it is a single $$[3\times 3]$$ matrix, that is
  a sum of all symmetry allowed elemenets multiplied by the given
  prefactors.
  
### See Also
  
[spinw.setmatrix](spinw_setmatrix)
 

{% include links.html %}
