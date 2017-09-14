---
{title: spinw.genlattice( ), summary: generates crystal lattice from given parameters,
  keywords: sample, sidebar: sw_sidebar, permalink: spinw_genlattice.html, folder: spinw,
  mathjax: 'true'}

---
 
{R} = GENLATTICE(obj, 'option1', value1, 'option2', value2...)
 
Input:
 
obj       spinw class object.
 
Options:
 
angled    Alpha, beta, gamma angles in degree, dimensions are [1 3].
angle     Alpha, beta, gamma angles in radian, dimensions are [1 3].
lat_const a, b, c lattice parameters, dimensions are [1 3].
spgr      Space group index, or space group label (string), or space group
          operators in a matrix with dimensions [3 4 nOp].
label     Optional label for the space group if the generators are given
          in the 'spgr' option.
bv        Basis vectors given in a matrix with dimensions of [3 3], where
          each column defines a basis vector.
origin    Origin for the space group operators. Default is [0 0 0].
perm      Permutation of the abc axes of the space group operators.
nformula  Gives the number of formula units in the unit cell. It is used
          to normalize cross section in absolute units. Default is 0,
          when cross section is normalized per unit cell.
 
Output:
 
R         Rotation matrix that brings the input basis vector to the SpinW
          compatible form:
                  BVspinw = R*BV
 
Alternatively the lattice parameters can be given directly when the spinw
object is created using: spinw(inpStr), where struct contains the fields
with initial parameters, e.g.:
  inpStr.lattice.lat_const = [3 3 4];
 
The sym option points to the appropriate line in the symmetry.dat file,
where every line defines a space group by its generators. The first 230
lines contains all crystallographic space groups with standard setting
as it is in the International Tables of Crystallography. Additional lines
can be added to the symmetry.dat file using the sw_addsym() function.
Every line in the symmetry.dat file can be referenced by either its line
index or by its label (string).
 
If the sym option is 0, no symmetry will be used. The spinw.gencoupling()
function will determine the equivalent bonds based on bond length.
 
Output:
 
The obj.lattice field will be changed based on the input, the lattice
constants stored directly and the optional space group string is
converted to the integer type index.
 
Example:
 
...
crystal.genlattice('lat_const',[3 3 4],'angled',[90 90 120],'spgr','P 6')
crystal.genlattice('lat_const',[3 3 4],'angled',[90 90 120],'spgr',168)
crystal.genlattice('lat_const',[3 3 4],'angled',[90 90 120],'spgr','-y,x-y,z; -x,-y,z','label','R -3 m')
 
The three lines are equivalent, both will create hexagonal lattice, with
'P 6' space group.
 
See also SPINW, SWSYM.ADD, SWSYM.OPERATOR, SPINW.GENCOUPLING.
 

