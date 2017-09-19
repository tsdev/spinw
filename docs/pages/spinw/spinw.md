---
{title: spinw class, link: spinw class, summary: class to store and solve magnetic
    Hamiltonians, keywords: sample, sidebar: sw_sidebar, permalink: spinw.html, folder: spinw,
  mathjax: 'true'}

---
 
* * *
`obj = spinw()`
* * *
 
constructs a new spinw class object, with default parameters.
 
* * *
`obj = spinw(obj)`
* * *
 
constructs new spinw class object. If `obj` is spinw class, it only
checks the integrity of its internal data structure. If `obj` is
struct type, it creates new spinw object and checks data integrity.
 
* * *
`obj = spinw(source)`
* * *
 
construct new spinw class object, where `source` is either a file
path pointing to a local cif or fst file or a link to an online file.
 
* * *
`obj = spinw(figure_handle)`
* * *
 
copy the spinw object stored in a previous structural3D plot figure.
 
The data structure behind the spinw object can be accessed by using
`struct(obj)`. All fields of the struct type data behind the spinw
object are accessible through the main field names of the `obj`
object. For example the lattice parameters can be extracted using:
```matlab
  abc = obj.unit_cell.lat_const
```
 
spinw is a handle class, that means that only the handle of the
object is copied in a `swobj1 = swobj2` command. To create a copy
(clone) of an spinw object use
```matlab
   swobj1 = swobj2.copy
```
 
### Methods
 
#### Lattice operations
 
  spinw.genlattice
  spinw.basisvector
  spinw.rl
  spinw.nosym
  spinw.newcell
* [spinw.addatom](spinw_addatom.html) adds new atom
  spinw.unitcell
  spinw.abc
  spinw.atom
  spinw.matom
  spinw.natom
* [spinw.formula](spinw_formula.html) returns chemical formula, mass, volume, etc.
  spinw.disp
  spinw.symmetry
    
#### Plotting
 
  spinw.plot
 
#### Crystallographic twin operations
 
  spinw.addtwin
  spinw.twinq
  spinw.notwin
  spinw.ntwin
 
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
  spinw.nmat
    
#### Spin Hamiltonian generations
 
  spinw.quickham
  spinw.gencoupling
  spinw.addcoupling
  spinw.couplingtable
  spinw.addaniso
  spinw.addg
  spinw.field
  spinw.nbond
  spinw.temperature
  spinw.intmatrix
  spinw.symop
  spinw.setunit
    
#### Calculators
 
* [spinw.spinwave](spinw_spinwave.html) calculates spin correlation function using linear spin wave theory
* [spinw.powspec](spinw_powspec.html) calculates powder averaged spin wave spectra
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
 

