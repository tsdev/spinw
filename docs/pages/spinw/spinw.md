---
{title: spinw class, link: spinw class, summary: class to store and solve magnetic
    Hamiltonians, keywords: sample, sidebar: sw_sidebar, permalink: spinw, folder: spinw,
  mathjax: 'true'}

---
 
### Syntax
 
`obj = spinw`
 
`obj = spinw(obj)`
 
`obj = spinw(source)`
 
`obj = spinw(figure_handle)`
 
### Description
 
`obj = spinw` constructs an empty spinw class object.
 
`obj = spinw(obj)` constructs a spinw class object from the
parameters defined in `obj`. If `obj` is spinw class, it only checks
the integrity of its internal data structure. If `obj` is struct
type, it creates new spinw object and checks data integrity.
 
`obj = spinw(source)` construct new spinw class object, where
`source` is either a file path pointing to a local `cif` or `fst`
file or a link to an online file.
 
`obj = spinw(figure_handle)` copy the spinw object stored in a
previous structural 3D plot figure, referenced by `figure_handle`.
 
 
The data structure within the spinw object can be accessed by using
[spinw.struct] method. All fields of the struct type data behind the
spinw object are accessible through the main field names of the `obj`
object. For example the lattice parameters can be accessed using:
 
```matlab
abc = obj.unit_cell.lat_const
```
 
spinw is a handle class, which means that only the handle of the
object is copied in an assinment command `swobj1 = swobj2`. To create
a copy (clone) of an spinw object use:
 
```matlab
swobj1 = swobj2.copy
```
 
### Properties
 
The data within the `spinw` object is organized into a tree structure
with the main groups and the type of data they store are the
following:
 
* [spinw.lattice] unit cell parameters
* [spinw.unit_cell] atoms in the crystallographic unit cell
* [spinw.twin] crystal twin parameters
* [spinw.matrix] 3x3 matrices for using them in the Hailtonian
* [spinw.single_ion] single ion terms of the Hamiltonian
* [spinw.coupling] list of bonds
* [spinw.mag_str] magnetic structure
* [spinw.unit] physical units for the Hamiltonian
* [spinw.cache] temporary values
 
### Methods
 
Methods are the different commands that require a `spinw` object as a
first input, thus they can be called as `method1(obj,...)`,
alternatively the equivalent command is `obj.method1(...)`. The list
of public methods is below.
 
#### Lattice operations
 
  spinw.genlattice
  spinw.basisvector
  spinw.rl
  spinw.nosym
  spinw.newcell
  spinw.addatom
  spinw.unitcell
  spinw.abc
  spinw.atom
  spinw.matom
  spinw.natom
  spinw.formula
  spinw.disp
  spinw.symmetry
    
#### Plotting
 
  spinw.plot
 
#### Crystallographic twin operations
 
  spinw.addtwin
  spinw.twinq
  spinw.notwin
 
#### Magnetic structure operations
 
  spinw.genmagstr
  spinw.magstr
  spinw.magtable
  spinw.nmagext
  spinw.optmagstr
  spinw.optmagk
  spinw.optmagsteep
  spinw.anneal
  spinw.annealloop
  spinw.structfact
    
#### Matrix operations
 
  spinw.addmatrix
  spinw.getmatrix
  spinw.setmatrix
    
#### Spin Hamiltonian generations
 
  spinw.quickham
  spinw.gencoupling
  spinw.addcoupling
  spinw.couplingtable
  spinw.addaniso
  spinw.addg
  spinw.field
  spinw.temperature
  spinw.intmatrix
  spinw.symop
  spinw.setunit
    
#### Solvers
 
  spinw.spinwave
  spinw.powspec
  spinw.energy
  spinw.moment
  spinw.spinwavesym
  spinw.symbolic
  spinw.meanfield
  spinw.fourier
  spinw.fouriersym
 
#### Fitting spin wave spectrum
 
  spinw.fitspec
  spinw.matparser
  spinw.horace
    
#### Miscellaneous
 
  spinw.copy
  spinw.export
  spinw.fileid
  spinw.table
  spinw.validate
  spinw.version
  spinw.struct
  spinw.clearcache
  spinw.spinw
 
### See also
 
[spinw.copy], [spinw.struct], [Comparing handle and value classes](https://www.google.ch/url?sa=t&rct=j&q=&esrc=s&source=web&cd=3&ved=0ahUKEwjCvbbctqTWAhVBblAKHQxnAnIQFggyMAI&url=https%3A%2F%2Fwww.mathworks.com%2Fhelp%2Fmatlab%2Fmatlab_oop%2Fcomparing-handle-and-value-classes.html&usg=AFQjCNFoN4qQdn6rPXKWkQ7aoog9G-nHgA)
 

{% include links.html %}
