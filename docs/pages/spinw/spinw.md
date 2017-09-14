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
 
* [spinw.genlattice](spinw_genlattice.html) generates crystal lattice
* [spinw.basisvector](spinw_basisvector.html) generates basis vectors of the crystal lattice
* [spinw.rl](spinw_rl.html) generates reciprocal lattice basis vectors of the crystal lattice
* [spinw.nosym](spinw_nosym.html) removes the space group symmetry
* [spinw.newcell](spinw_newcell.html) changes lattice vectors while keeping atoms
* [spinw.addatom](spinw_addatom.html) adds new atom to an spinw object
* [spinw.unitcell](spinw_unitcell.html) returns information on atoms in the crystallographic unit cell
* [spinw.abc](spinw_abc.html) returns lattice parameters and angles
* [spinw.atom](spinw_atom.html) generates all atomic positions in the unit cell
* [spinw.matom](spinw_matom.html) generates all magnetic atoms in the unit cell
* [spinw.natom](spinw_natom.html) gives the number of symmetry unrelated atoms in the unit cell
* [spinw.formula](spinw_formula.html) returns chemical formula, mass, volume, etc.
* [spinw.disp](spinw_disp.html) prints the spinw data structure in readable format onto the Command Window
* [spinw.symmetry](spinw_symmetry.html) true if space group is used to generate couplings
    
#### Plotting
 
* [spinw.plot](spinw_plot.html) plots crystal structure, magnetic structure, anisotropy and bonds
 
#### Crystallographic twin operations
 
* [spinw.addtwin](spinw_addtwin.html) adds new twins to an spinw object
* [spinw.twinq](spinw_twinq.html) calculates equivalent Q point in twins
* [spinw.notwin](spinw_notwin.html) removes any twin added to the spinw object
* [spinw.ntwin](spinw_ntwin.html) gives the number of twins
 
#### Magnetic structure operations
 
* [spinw.genmagstr](spinw_genmagstr.html) generates magnetic structure
* [spinw.magstr](spinw_magstr.html) generates magnetic structure for the rotating frame
* [spinw.magtable](spinw_magtable.html) creates tabulated list of all magnetic moments stored in obj
* [spinw.nmagext](spinw_nmagext.html) gives the number of magnetic atoms in the magnetic supercell
* [spinw.optmagstr](spinw_optmagstr.html) optimises magnetic structure by minimizing the energy using non-linear optimization algorithms
* [spinw.optmagk](spinw_optmagk.html) determines the magnetic propagation vector
* [spinw.optmagsteep](spinw_optmagsteep.html) optimise magnetic structure using the method of steepest descent
* [spinw.anneal](spinw_anneal.html) performs simulated annealing on the magnetic structure
* [spinw.annealloop](spinw_annealloop.html) performs simulated annealing on the magnetic structure and measurements
* [spinw.structfact](spinw_structfact.html) calculates magnetic and nuclear structure factor
    
#### Matrix operations
 
* [spinw.addmatrix](spinw_addmatrix.html) adds new matrix that can be assigned to spins in the Hamiltonian
* [spinw.getmatrix](spinw_getmatrix.html) gives the symmetry allowed matrices for a given coupling or anisotropy
* [spinw.setmatrix](spinw_setmatrix.html) changes the selected matrix inside the spinw object
* [spinw.nmat](spinw_nmat.html) gives the number of matrices defined in an spinw object
    
#### Spin Hamiltonian generations
 
* [spinw.quickham](spinw_quickham.html) creates magnetic Hamiltonian with a single command
* [spinw.gencoupling](spinw_gencoupling.html) generates the COUPLING property of spinw object
* [spinw.addcoupling](spinw_addcoupling.html) assigns a predefined matrix as exchange coupling on selected bonds
* [spinw.couplingtable](spinw_couplingtable.html) creates tabulated list of all bonds as stored
* [spinw.addaniso](spinw_addaniso.html) assigns anisotropy matrices to magnetic ions
* [spinw.addg](spinw_addg.html) assigns g-tensor to magnetic ions
* [spinw.field](spinw_field.html) get/set magnetic field value
* [spinw.nbond](spinw_nbond.html) gives the number of bonds defined in the spinw object
* [spinw.temperature](spinw_temperature.html) get/set stored temperature value
* [spinw.intmatrix](spinw_intmatrix.html) creates the interactions matrices (connectors and values)
* [spinw.symop](spinw_symop.html) generates the symmetry operators on bonds and magnetic atoms
* [spinw.setunit](spinw_setunit.html) sets the physical units
    
#### Calculators
 
* [spinw.spinwave](spinw_spinwave.html) calculates spin correlation function using linear spin wave theory
* [spinw.powspec](spinw_powspec.html) calculates powder averaged spin wave spectra
* [spinw.energy](spinw_energy.html) calculates the ground state energy per spin
* [spinw.moment](spinw_moment.html) calculates the size of the reduced moment due to quantum and thermal fluctuations
* [spinw.spinwavesym](spinw_spinwavesym.html) calculates symbolic spin wave dispersion
* [spinw.symbolic](spinw_symbolic.html) change between symbolic/numeric calculation
* [spinw.meanfield](spinw_meanfield.html) mean field calculation of the wave vector dependent susceptibility
* [spinw.fourier](spinw_fourier.html) calculates the Fourier transformation of the Hamiltonian
* [spinw.fouriersym](spinw_fouriersym.html) calculates the Fourier transformation of a symbolic Hamiltonian
 
#### Fitting spin wave spectrum
 
* [spinw.fitspec](spinw_fitspec.html) fits spin wave spectra to experimental spectral data
* [spinw.matparser](spinw_matparser.html) assigns new values to existing matrices
* [spinw.horace](spinw_horace.html) calculates spin wave dispersion/correlation functions to be called from Horace
    
#### Miscellaneous
 
* [spinw.copy](spinw_copy.html) clones spinw object with all data
* [spinw.export](spinw_export.html) export data from spinw object into different file formats
* [spinw.fileid](spinw_fileid.html) determines where the text out is written
* [spinw.table](spinw_table.html) outputs easy to read tables of internal data
* [spinw.validate](spinw_validate.html) validates spinw object properties
* [spinw.version](spinw_version.html) returns the version of SpinW used to create the object
* [spinw.struct](spinw_struct.html) extracts all public properties of spinw object into a struct
* [spinw.clearcache](spinw_clearcache.html) clears the cache
* [spinw.spinw](spinw_spinw.html) SPINW constructor
 
### Notes
 
Tutorials and documentation can be found at [psi.ch/spinw](https://psi.ch/spinw)
 
Forum for questions on [Google Groups](https://groups.google.com/forum/#!forum/spinwforum)
 
Lates version and bug reports/feature requests can be submitted on [GitHub](https://github.com/tsdev/spinw)
 
### See also
 
[spinw.copy](spinw_copy.html), [spinw.struct](spinw_struct.html), [Comparing handle and value classes](https://www.google.ch/url?sa=t&rct=j&q=&esrc=s&source=web&cd=3&ved=0ahUKEwjCvbbctqTWAhVBblAKHQxnAnIQFggyMAI&url=https%3A%2F%2Fwww.mathworks.com%2Fhelp%2Fmatlab%2Fmatlab_oop%2Fcomparing-handle-and-value-classes.html&usg=AFQjCNFoN4qQdn6rPXKWkQ7aoog9G-nHgA)
 

