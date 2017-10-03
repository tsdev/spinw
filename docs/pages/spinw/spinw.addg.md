---
{title: spinw.addg method, link: spinw.addg, summary: assigns g-tensor to magnetic
    atoms, keywords: sample, sidebar: sw_sidebar, permalink: spinw_addg, folder: spinw,
  mathjax: 'true'}

---
  
### Syntax
  
`addg(obj, matrixIdx, {atomTypeIdx}, {atomIdx})`
  
### Description
  
`addg(obj, matrixIdx, {atomTypeIdx}, {atomIdx})` assigns the
$$[3\times 3]$$ matrix selected by `matrixIdx` (using either the matrix
label or matrix index) to the magnetic sites selected by `atomTypeIdx`
that can contain a name of an atom or its atom index (see [spinw.atom](spinw_atom)).
If `atomTypeIdx` is not defined, g-tensor will be assigned to all
magnetic atoms.
  
### Examples
  
The following example will add the $$g_1$$ diagonal matrix to all magnetic
atoms as anisotropic g-tensor:
  
```matlab
cryst = spinw
cryst.genlattice('lat_const',[4 4 3],'spgr','P 4')
cryst.addatom('r',[1/4 1/4 1/2],'S',1)
cryst.addmatrix('label','g_1','value',diag([-0.1 0 0]))
cryst.gencoupling
cryst.addg('g_1')
cryst.plot('ionMode','g')
```
 
{% include image.html file="generated/spinw_1.png" alt="cryst.plot('ionMode','g')" %}
  
### Input Arguments
  
`matrixIdx`
: Either an integer, that selects the matrix
  `obj.matrix.mat(:,:,matrixIdx)`, or a string identical to one
  of the previously defined matrix labels, stored in
  `obj.matrix.label`. Maximum value is $$n_{mat}$$.
  
`atomTypeIdx`
: String or cell of strings that select magnetic atoms by
  their label. Also can be a vector that contains integers, the index of
  the magnetic atoms in [spinw.unit_cell](spinw_unit_cell), this will assign the given
  g-tensor to all symmetry equivalent atoms. Maximum value is $$n_{atom}$$.
  If `atomTypeIdx` is not defined, the given g-tensor will be assigned to
  all magnetic atoms. Optional.
 
`atomIdx`
: A vector that contains indices selecting some of the
  symmetry equivalent atoms. Maximum value is the number of symmetry
  equivalent atoms corresponding to `atomTypeIdx`. If the crystal
  symmetry is higher than $$P0$$, `atomIdx` is not allowed, since the
  g-tensor for equivalent atoms will be calculated using the symmetry
  operators of the space group. Optional.
  
### Output Arguments
  
The function adds extra entries to the `obj.single_ion.g` matrix.
  
### See Also
  
[spinw](spinw) \| [spinw.addcoupling](spinw_addcoupling) \| [spinw.addaniso](spinw_addaniso) \| [spinw.addmatrix](spinw_addmatrix)
 

{% include links.html %}
