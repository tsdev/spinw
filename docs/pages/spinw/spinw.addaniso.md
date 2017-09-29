---
{title: spinw.addaniso method, link: spinw.addaniso, summary: assigns anisotropy to
    magnetic sites, keywords: sample, sidebar: sw_sidebar, permalink: spinw_addaniso.html,
  folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`addaniso(obj, matrixIdx, {atomTypeIdx}, {atomIdx})`
  
### Description
  
`addaniso(obj, matrixIdx, {atomTypeIdx}, {atomIdx})` assigns the
$$[3\times 3]$$ matrix selected by `matrixIdx` (using either the matrix
label or matrix index) to the magnetic sites selected by `atomTypeIdx`
that can contain a name of an atom or its atom index (see [spinw.atom](spinw_atom.html)).
If `atomTypeIdx` is not defined, anisotropy will be assigned to all
magnetic atoms.
  
### Examples
  
To add the $$A_1$$ diagonal matrix to all magnetic atoms as
anisotropy (easy *XY* plane anisotropy) use the following code:
 
```matlab
cryst = spinw
cryst.genlattice('lat_const',[4 4 3],'spgr','P 4')
cryst.addatom('r',[1/4 1/4 1/2],'S',1)
cryst.addmatrix('label','A1','value',diag([-0.1 0 0]))
cryst.gencoupling
cryst.addaniso('A1')
plot(cryst)
```
 
{% include image.html file="generated/spinw_addaniso_1.png" alt="plot(cryst)" %}
  
### Input arguments
 
`obj`
: [spinw](spinw.html) object.
 
### Name-Value Pair Arguments
  
`matrixIdx`
: Either an integer, that selects the matrix according to
  `obj.matrix.mat(:,:,matrixIdx)`, or a string identical to one
  of the previously defined matrix labels, stored in
  `obj.matrix.label`.
  
`atomTypeIdx`
: String or cell of strings that select magnetic atoms by
  their label. Also can be a vector that contains integers, the index of
  the magnetic atoms in `obj.unit_cell`, with all symmetry equivalent
  atoms. Maximum value is $$n_{atom}$$, if undefined anisotropy is assigned to
  all magnetic atoms. Optional.
 
`atomIdx`
: A vector that contains indices selecting some of the
  symmetry equivalent atoms. Maximum value is the number of symmetry
  equivalent atoms generated corresponding to `atomTypeIdx` site. If
  crystal symmetry is not 0, `atomIdx` is not allowed, since the
  anisotropy matrix for equivalent atoms will be calculated using the
  symmetry operators of the space group. Optional.
  
### Output Arguments
  
The function adds extra entries in the `obj.single_ion.aniso` field of the
obj [spinw](spinw.html) object.
  
### See Also
  
[spinw](spinw.html), [spinw.single_ion](spinw_single_ion.html), [spinw.addcoupling](spinw_addcoupling.html), [spinw.addg](spinw_addg.html) and [spinw.addmatrix](spinw_addmatrix.html)
 

{% include links.html %}
