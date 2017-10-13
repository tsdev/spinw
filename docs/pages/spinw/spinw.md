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
[spinw.struct](spinw_struct) method. All fields of the struct type data behind the
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
 
* [spinw.lattice](spinw_lattice) unit cell parameters
* [spinw.unit_cell](spinw_unit_cell) atoms in the crystallographic unit cell
* [spinw.twin](spinw_twin) crystal twin parameters
* [spinw.matrix](spinw_matrix) 3x3 matrices for using them in the Hailtonian
* [spinw.single_ion](spinw_single_ion) single ion terms of the Hamiltonian
* [spinw.coupling](spinw_coupling) list of bonds
* [spinw.mag_str](spinw_mag_str) magnetic structure
* [spinw.unit](spinw_unit) physical units for the Hamiltonian
* [spinw.cache](spinw_cache) temporary values
 
### Methods
 
Methods are the different commands that require a `spinw` object as a
first input, thus they can be called as `method1(obj,...)`,
alternatively the equivalent command is `obj.method1(...)`. The list
of public methods is below.
 
#### Lattice operations
 
* [spinw.genlattice](spinw_genlattice) generates crystal lattice
* [spinw.basisvector](spinw_basisvector) generates lattice vectors
* [spinw.rl](spinw_rl) generates reciprocal lattice vectors
* [spinw.nosym](spinw_nosym) reduces symmetry to P0
* [spinw.newcell](spinw_newcell) transforms lattice
* [spinw.addatom](spinw_addatom) adds new atom
* [spinw.unitcell](spinw_unitcell) returns unit cell data
* [spinw.abc](spinw_abc) returns lattice parameters and angles
* [spinw.atom](spinw_atom) generates symmetry equivalent atomic positions
* [spinw.matom](spinw_matom) generates magnetic lattice
* [spinw.natom](spinw_natom) number of symmetry unrelated atoms
* [spinw.formula](spinw_formula) returns basic physical properties
* [spinw.disp](spinw_disp) prints information
* [spinw.symmetry](spinw_symmetry) returns whether symmetry is defined
    
#### Plotting
 
* [spinw.plot](spinw_plot) plots 3D model
 
#### Crystallographic twin operations
 
* [spinw.addtwin](spinw_addtwin) adds crystallographic twins
* [spinw.twinq](spinw_twinq) calculates equivalent Q point in twins
* [spinw.notwin](spinw_notwin) removes all twins
 
#### Magnetic structure operations
 
* [spinw.genmagstr](spinw_genmagstr) generates magnetic structure
* [spinw.magstr](spinw_magstr) returns single-k magnetic structure representation
* [spinw.magtable](spinw_magtable) creates tabulated list of all magnetic moments stored in obj
* [spinw.nmagext](spinw_nmagext) number of magnetic sites
* [spinw.optmagstr](spinw_optmagstr) general magnetic structure optimizer
* [spinw.optmagk](spinw_optmagk) determines the magnetic propagation vector
* [spinw.optmagsteep](spinw_optmagsteep) quench optimization of magnetic structure
* [spinw.anneal](spinw_anneal) performs simulated annealing of spins
* [spinw.annealloop](spinw_annealloop) parameter sweep for simulated annealing
* [spinw.structfact](spinw_structfact) calculates magnetic and nuclear structure factor
    
#### Matrix operations
 
* [spinw.addmatrix](spinw_addmatrix) adds new [3x3] matrix
* [spinw.getmatrix](spinw_getmatrix) determines the symmetry allowed tensor elements
* [spinw.setmatrix](spinw_setmatrix) sets exchange tensor values
    
#### Spin Hamiltonian generations
 
* [spinw.quickham](spinw_quickham) quickly generate magnetic Hamiltonian
* [spinw.gencoupling](spinw_gencoupling) generates bond list
* [spinw.addcoupling](spinw_addcoupling) assigns an exchange matrix to a bond
  spinw.couplingtable
* [spinw.addaniso](spinw_addaniso) assigns anisotropy to magnetic sites
* [spinw.addg](spinw_addg) assigns g-tensor to magnetic atoms
* [spinw.field](spinw_field) get/set magnetic field value
* [spinw.temperature](spinw_temperature) get/set temperature
* [spinw.intmatrix](spinw_intmatrix) generates interaction matrix
* [spinw.symop](spinw_symop) generates the bond symmetry operators
* [spinw.setunit](spinw_setunit) sets the physical units
    
#### Solvers
 
* [spinw.spinwave](spinw_spinwave) calculates spin correlation function using linear spin wave theory
* [spinw.powspec](spinw_powspec) calculates powder averaged spin wave spectra
* [spinw.energy](spinw_energy) calculates the ground state energy
* [spinw.moment](spinw_moment) calculates quantum correction on ordered moment
* [spinw.spinwavesym](spinw_spinwavesym) calculates symbolic spin wave dispersion
* [spinw.symbolic](spinw_symbolic) switches between symbolic/numeric mode
  spinw.meanfield
* [spinw.fourier](spinw_fourier) calculates the Fourier transformation of the Hamiltonian
* [spinw.fouriersym](spinw_fouriersym) calculates the Fourier transformation of the symbolic Hamiltonian
 
#### Fitting spin wave spectrum
 
* [spinw.fitspec](spinw_fitspec) fits experimental spin wave data
* [spinw.matparser](spinw_matparser) parses parameter vector into matrices
* [spinw.horace](spinw_horace) spin wave calculator with interface to Horace
    
#### Miscellaneous
 
* [spinw.copy](spinw_copy) clones spinw object
* [spinw.export](spinw_export) export data into file
* [spinw.fileid](spinw_fileid) determines file object for text output
* [spinw.table](spinw_table) outputs easy to read tables of internal data
* [spinw.validate](spinw_validate) validates spinw object properties
* [spinw.version](spinw_version) returns the version of SpinW
* [spinw.struct](spinw_struct) converts properties into struct
* [spinw.clearcache](spinw_clearcache) clears the cache
* [spinw.spinw](spinw_spinw) Operators and special characters.
 
### See also
 
[spinw.copy](spinw_copy), [spinw.struct](spinw_struct), [Comparing handle and value classes](https://www.google.ch/url?sa=t&rct=j&q=&esrc=s&source=web&cd=3&ved=0ahUKEwjCvbbctqTWAhVBblAKHQxnAnIQFggyMAI&url=https%3A%2F%2Fwww.mathworks.com%2Fhelp%2Fmatlab%2Fmatlab_oop%2Fcomparing-handle-and-value-classes.html&usg=AFQjCNFoN4qQdn6rPXKWkQ7aoog9G-nHgA)
 

{% include links.html %}
