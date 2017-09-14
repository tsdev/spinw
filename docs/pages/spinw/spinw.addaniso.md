---
{title: spinw.addaniso( ), summary: assigns anisotropy matrices to magnetic ions,
  keywords: sample, sidebar: sw_sidebar, permalink: spinw_addaniso.html, folder: spinw,
  mathjax: 'true'}

---
 
ADDANISO(obj, matrixIdx, {atomTypeIdx}, {atomIdx})
 
Input:
 
matrixIdx     Either an integer, that selects the matrix
              obj.matrix.mat(:,:,matrixIdx), or a string identical to one
              of the previously defined matrix labels, stored in
              obj.matrix.label. Maximum value is nJ.
atomTypeIdx   String or cell of strings that select magnetic atoms by
              their label. Also can be a vector that contains integers,
              the index of the magnetic atoms in obj.unit_cell, with all
              symmetry equivalent atoms. Maximum value is nAtom, if
              undefined anisotropy is assigned to all magnetic atoms.
              Optional.
 atomIdx      A vector that contains indices selecting some of the
              symmetry equivalent atoms. Maximum value is the number of
              symmetry equivalent atoms generated. If crystal symmetry is
              not 0, atomIdx is not allowed, since the anisotropy matrix
              for equivalent atoms will be calculated using the symmetry
              operators of the space group. Optional.
 
Output:
 
The function adds extra entries in the 'single_ion.aniso' field of the
obj spinw object.
 
Example:
 
...
cryst.addmatrix('label','A1','value',diag([-0.1 -0.1 0]))
cryst.gencoupling
cryst.addaniso('A1')
 
This will add the 'A1' diagonal matrix to all magnetic atoms as
anisotropy (easy XY plane anisotropy).
 
See also SPINW, SPINW.ADDCOUPLING, SPINW.ADDG, SPINW.ADDMATRIX.
 

