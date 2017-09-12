---
{title: '@spinw.newcell( )', summary: changes lattice vectors while keeping atoms,
  keywords: sample, sidebar: sw_sidebar, permalink: spinw_newcell.html, folder: '@spinw',
  mathjax: 'true'}

---
changes lattice vectors while keeping atoms
 
{T} = NEWCELL(obj, 'option1', value1, ...)
 
The function defines new unit cell using the 3 vectors contained in
bvect. The three vectors in lattice units define a parallelepiped. This
will be the new unit cell. The atoms from the original unit cell will
fill the new unit cell. Also the magnetic structure and bond and single
ion property definitions will be erased from the structure. The new cell
will naturally have different reciprocal lattice, however the original
reciprocal lattice units will be retained automatically, to use the new
reciprocal lattice, set the 'keepq' option to false. In the default case
the spinw.spinwave() function will calculate spin wave dispersion at
reciprocal lattice points of the original lattice. The transformation
between the two lattices is stored in spinw.unit.qmat.
 
Input:
 
obj       spinw class object.
 
Options:
 
bvect     Defines the new lattice vectors in the original lattice
          coordinate system. Cell with the following elements
          {v1 v2 v3} or a 3x3 matrix with v1, v2 and v3 as column
          vectors: [v1 v2 v3].
bshift    Row vector defines a shift of the position of the unit cell.
          Optional.
keepq     If true, the reciprocal lattice units of the new model will be
          the same as in the old model. This is achieved by storing the
          transformation matrix between the new and the old coordinate
          system in spinw.unit.qmat. Default is true.
 
Output:
 
T     is a transformation matrix that converts Q points (in reciprocal
      lattice units) from the old reciprocal lattice to the new
      reciprocal lattice as follows:
          Qrlu_new = T * Qrlu_old,
      where the dimensions of the Q vectors are [1 3].
 
Example:
 
In this example we generate the triangular lattice antiferromagnet and
convert the hexagonal cell to orthorhombic. This doubles the number of
magnetic atoms in the cell and changes the reciprocal lattice. However we
use 'keepq' to able to index the reciprocal lattice of the orthorhombic
cell with the reciprocal lattice of the original hexagonal cell. To show
that the two models are equivalent, we calculate the spin wave spectrum
on both model using the same rlu. On the orthorhombic cell, the q value
will be converted automatically and the calculated spectrum will be the
same for both cases.
 
tri = sw_model('triAF',1);
tri_orth = copy(tri);
tri_orth.newcell('bvect',{[1 0 0] [1 2 0] [0 0 1]},'keepq',true);
tri_orth.gencoupling
tri_orth.addcoupling('bond',1,'mat','J1')
newk = ((tri_orth.unit.qmat)*tri.magstr.k')';
tri_orth.genmagstr('mode','helical','k',newk,'S',[1 0 0]')
plot(tri_orth)
  
figure
subplot(2,1,1)
sw_plotspec(sw_egrid(tri.spinwave({[0 0 0] [1 1 0] 501})),'mode','color','dE',0.2)
subplot(2,1,2)
sw_plotspec(sw_egrid(tri_orth.spinwave({[0 0 0] [1 1 0] 501})),'mode','color','dE',0.2)
 
 
See also SPINW.GENLATTICE, SPINW.GENCOUPLING, SPINW.NOSYM.
 
