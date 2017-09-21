---
{title: spinw.genlattice method, link: spinw.genlattice, summary: generates crystal
    lattice, keywords: sample, sidebar: sw_sidebar, permalink: spinw_genlattice.html,
  folder: spinw, mathjax: 'true'}

---
  
`{R} = genlattice(obj, 'option1', value1, ...)`
* * *
 
The function generates all necessary parameters to define a lattice with
symmetry and store it in `obj.lattice`.
 
Alternatively the lattice parameters can be given directly when the spinw
object is created using `spinw(inpStr)` command, where struct contains
the fields with initial parameters, e.g.:
```matlab
inpStr.lattice.lat_const = [3 3 4];
```
 
### Input
 
`obj`
: [spinw](spinw.html) object.
  
### Options
  
`angled`
: `[α, β, γ]` angles in °, dimensions are $$[1\times 3]$$.
  
`angle`
: `[α, β, γ]` angles in radian, dimensions are $$[1\times 3]$$.
  
`lat_const`
: `[a, b, c]` lattice parameters in units defined in [spinw.unit](spinw_unit.html) (with Å
  being the default), dimensions are $$[1\times 3]$$.
  
`spgr`
: Defines the space group. Can have the following values:
 
  * **space group label** string, name of the space group, can be any
    label defined in the `symmetry.dat` file.
  * **space group index** line number in the `symmetry.dat` file.
  * **space group operators** matrix with dimensions 
    $$[3\times 4\times n_{op}]$$.
    
  The `symmetry.dat` file stores definition of the 230 space groups in
  standard settings as it is in the [International Tables of Crystallography](http://it.iucr.org/A/).
  Additional lines can be added to the `symmetry.dat` file using the
  [swsym.add](swsym_add.html) function which later can be used in the `spgr` option.
  
  If the `spgr` option is 0, no symmetry will be used. The
  [spinw.gencoupling](spinw_gencoupling.html) function will determine the equivalent bonds based on
  bond length.
  
`label`
: Optional label for the space group if the generators are given in the
  `spgr` option.
`bv`
: Basis vectors given in a matrix with dimensions of $$[3\times 3]$$, where
  each column defines a basis vector.
  
`origin`
: Origin for the space group operators, default value is `[0 0 0]`.
  
`perm`
: Permutation of the abc axes of the space group operators.
  
`nformula`
: Gives the number of formula units in the unit cell. It is used
  to normalize cross section in absolute units. Default value is 0, when
  cross section is normalized per unit cell.
  
### Output
  
`R`
: Rotation matrix that brings the input basis vector to the SpinW
  compatible form:
  ```matlab
  BVspinw = R*BV
  ```
  
The result of the `spinw.genlattice` function is that `obj.lattice` field
will be changed based on the input, the lattice parameters are stored
directly and the optional space group string is converted into space
group operator matrices.
 
### Example
 
```matlab
crystal.genlattice('lat_const',[3 3 4],'angled',[90 90 120],'spgr','P 6')
crystal.genlattice('lat_const',[3 3 4],'angled',[90 90 120],'spgr',168)
crystal.genlattice('lat_const',[3 3 4],'angled',[90 90 120],'spgr','-y,x-y,z; -x,-y,z','label','R -3 m')
```
 
The three lines are equivalent, both will create hexagonal lattice, with
$$P6$$ space group.
 
### See also
 
[spinw](spinw.html), [swsym.add](swsym_add.html), [swsym.operator](swsym_operator.html), [spinw.gencoupling](spinw_gencoupling.html)
 

{% include links.html %}
