---
{title: spinw.newcell method, link: spinw.newcell, summary: transforms lattice, keywords: sample,
  sidebar: sw_sidebar, permalink: spinw_newcell, folder: spinw, mathjax: 'true'}

---
  
### Syntax
  
`{t} = newcell(obj,Name,Value)`
  
### Description
  
`{t} = newcell(obj,Name,Value)` redefines the unit cell using new basis
vectors. The input three basis vectors are in lattice units of the
original cell and define a parallelepiped. The atoms from the original
unit cell will fill the new unit cell and if the two cells are compatible
the structure won't change. The magnetic structure, bonds and single ion
property definitions will be erased. The new cell will have different
reciprocal lattice, however the original reciprocal lattice units will be
retained automatically. To use the new reciprocal lattice, set the
`'keepq'` option to `false`. In the default case the [spinw.spinwave](spinw_spinwave)
function will calculate spin wave dispersion at reciprocal lattice points
of the original lattice. The transformation between the two lattices is
stored in `spinw.unit.qmat`.
  
### Examples
  
In this example we generate the triangular lattice antiferromagnet and
convert the hexagonal cell to orthorhombic. This doubles the number of
magnetic atoms in the cell and changes the reciprocal lattice. However we
set `'keepq'` parameter to `true` to able to index the reciprocal lattice
of the orthorhombic cell with the reciprocal lattice of the original
hexagonal cell. To show that the two models are equivalent, we calculate
the spin wave spectrum on both model using the same rlu. On the
orthorhombic cell, the $$Q$$ value will be converted automatically and the
calculated spectrum will be the same for both cases.
 
```matlab
tri = sw_model('triAF',1)
tri_orth = copy(tri)
tri_orth.newcell('bvect',{[1 0 0] [1 2 0] [0 0 1]},'keepq',true)
tri_orth.gencoupling
tri_orth.addcoupling('bond',1,'mat','J_1')
newk = ((tri_orth.unit.qmat)*tri.magstr.k')'
tri_orth.genmagstr('mode','helical','k',newk,'S',[1 0 0]')
plot(tri_orth)
```
 
{% include image.html file="generated/spinw_ne_1.png" alt="swplot.zoom(1.5)" %}
```matlab
subplot(2,1,1)
sw_plotspec(sw_egrid(tri.spinwave({[0 0 0] [1 1 0] 501})),'mode','color','dE',0.2)
subplot(2,1,2)
spec = tri_orth.spinwave({[0 0 0] [1 1 0] 501});
sw_plotspec(sw_egrid(tri_orth.spinwave({[0 0 0] [1 1 0] 501})),'mode','color','dE',0.2)
```
 
{% include image.html file="generated/spinw_ne_2.png" alt="sw_plotspec(sw_egrid(tri_orth.spinwave({[0 0 0] [1 1 0] 501})),'mode','color','dE',0.2)" %}
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
### Name-Value Pair Arguments
  
`'bvect'`
: Defines the new lattice vectors in the original lattice
  coordinate system. Cell with the following elements
  `{v1 v2 v3}` or a $$[3\times 3]$$ matrix with `v1`, `v2` and `v3` as column
  vectors: `[v1 v2 v3]`. Default value is `eye(3)` for indentity
  transformation.
  
`'bshift'`
: Row vector that defines a shift of the position of the unit cell.
  Default value is `[0 0 0]`.
  
`'keepq'`
: If true, the reciprocal lattice units of the new model will be
  the same as in the old model. This is achieved by storing the
  transformation matrix between the new and the old coordinate system in
  `spinw.unit.qmat` and applying it every time a reciprocal space
  definition is invoked, such as in [spinw.spinwave](spinw_spinwave). Default value is
  `false`.
  
### Output Arguments
  
`T`
: Transformation matrix that converts $$Q$$ points (in reciprocal
      lattice units) from the old reciprocal lattice to the new
      reciprocal lattice as follows:
  ```matlab
  Qrlu_new = T * Qrlu_old
  ```
  where the $$Q$$ vectors are row vectors with 3 elements.
  
### See Also
  
[spinw.genlattice](spinw_genlattice) \| [spinw.gencoupling](spinw_gencoupling) \| [spinw.nosym](spinw_nosym)
 
*[rlu]: reciprocal lattice unit
 

{% include links.html %}
