---
{title: Class spinw, summary: SPINW class defines data structure and methods to calculate
    spin wave, keywords: sample, sidebar: sw_sidebar, permalink: spinw.html, folder: spinw,
  mathjax: 'true'}

---
SPINW class defines data structure and methods to calculate spin wave
dispersion in magnetic crystals.
 
obj = SPINW()
 
constructs a new spinw class object, with default parameters.
 
obj = SPINW(obj)
 
constructs new spinw class object. If obj is spinw class, it only
checks its data integrity. If obj is struct type, it creates new
spinw object and checks data integrity.
 
obj = SPINW(file_path)
 
construct new spinw class object, where file_path contains a string
of a .cif/.fst file defining an input crystal and/or magnetic
structure.
 
obj = SPINW(figure_handle)
 
copy the SpinW object stored in a previous structural plot.
 
The data structure behind the spinw object can be accessed by using
STRUCT(obj). All fields of the struct type data behind the spinw
object are accessible through the main field names of the obj object.
For example the lattice parameters:
  abc = obj.unit_cell.lat_const;
 
spinw is a handle class, that means that only the handle of the
object is copied in a swobj1 = swobj2 command. To create a copy
(clone) of an spinw object use:
   swobj1 = swobj2.COPY;
 
See also: SPINW.COPY, SPINW.STRUCT,
<a href='/Applications/MATLAB_R2012b.app/help/matlab/matlab_oop/comparing-handle-and-value-classes.html'>Comparing handle and value classes</a>.
 
### Methods
 
#### Lattice operations

* [spinw.genlattice](/spinw_genlattice) defines lattice parameters, angles and symmetry
* [spinw.basisvector](/spinw_basisvector) returns the basis vectors in a matrix
* [spinw.rl](/spinw_rl) returns the reciprocal lattice vectors in a matrix
* [spinw.nosym](/spinw_nosym) reduces the lattice symmetry to P0
* [spinw.newcell](/spinw_newcell) generates an arbitrary new cell and reduces symmetry to P0
* [spinw.addatom](/spinw_addatom) adds a new atoms to the lattice
* [spinw.unitcell](/spinw_unitcell) selects a subset of the symmetry inequivalent atoms
* [spinw.abc](/spinw_abc) returns a vector [a b c alpha beta gamma]
* [spinw.atom](/spinw_atom) returns all symmetry generated atom
* [spinw.matom](/spinw_matom) returns all symmetry generated magnetic atoms
* [spinw.natom](/spinw_natom) returns the number of symmetry inequivalent atoms
* [spinw.formula](/spinw_formula) returns basic unit cell information
* [spinw.disp](/spinw_disp) return basic lattice information
* [spinw.symmetry](/spinw_symmetry) returns whether symmetry is on (>P0)
    
#### Plotting

* [spinw.plot](/spinw_plot) plots crystal, magnetic structure and spin Hamiltonian
 
#### Crystallographic twin operations

* [spinw.addtwin](/spinw_addtwin) adds crystallographic twin
* [spinw.twinq](/spinw_twinq) determines reciprocal lattice positions in twins
* [spinw.notwin](/spinw_notwin) removes the twins
* [spinw.ntwin](/spinw_ntwin) number of twins
 
#### Magnetic structure operations

* [spinw.genmagstr](/spinw_genmagstr) generate different types of magnetic structures
* [spinw.magstr](/spinw_magstr) returns the magnetic structure in a rotating frame representation
* [spinw.magtable](/spinw_magtable) returns a table of the magnetic structure
* [spinw.nmagext](/spinw_nmagext) returns the number of magnetic moments in the supercell
* [spinw.optmagstr](/spinw_optmagstr) optimizes magnetic structure with constraints
* [spinw.optmagk](/spinw_optmagk) optimizes the magnetic propagation vector
* [spinw.optmagsteep](/spinw_optmagsteep) optimizes the moment directions for fixed-k
* [spinw.anneal](/spinw_anneal) simulated annealing of magnetic structures
* [spinw.annealloop](/spinw_annealloop) simulated annealing and scanning a parameter in the Hamiltonian
* [spinw.structfact](/spinw_structfact) calculates magnetic structure factor (EXPERIMENTAL)
    
#### Matrix operations

* [spinw.addmatrix](/spinw_addmatrix) add a new matrix to the internal data
* [spinw.getmatrix](/spinw_getmatrix) determines the symmetry allowed elements of a matrix
* [spinw.setmatrix](/spinw_setmatrix) generates matrix values based on symmetry allowed elements
* [spinw.nmat](/spinw_nmat) number of matrices
    
#### Spin Hamiltonian generations

* [spinw.quickham](/spinw_quickham) quick generation of Heisenberg Hamiltonian
* [spinw.gencoupling](/spinw_gencoupling) generates the list of bonds
* [spinw.addcoupling](/spinw_addcoupling) assigns a matrix to a bond
* [spinw.couplingtable](/spinw_couplingtable) lists information on the bonds
* [spinw.addaniso](/spinw_addaniso) assigns a matrix to a magnetic ion as anisotropy
* [spinw.addg](/spinw_addg) assigns a matrix to a magnetic ion as g-tensor
* [spinw.field](/spinw_field) stores the magnetic field
* [spinw.nbond](/spinw_nbond) number of bonds
* [spinw.temperature](/spinw_temperature) temperature for thermal population calculation
* [spinw.intmatrix](/spinw_intmatrix) returns the spin Hamiltonian after symmetry applied (INTERNAL)
* [spinw.symop](/spinw_symop) generates the symmetry operators on bonds and magnetic atoms (INTERNAL)
* [spinw.setunit](/spinw_setunit) sets the physical units
    
#### Calculators

* [spinw.spinwave](/spinw_spinwave) linear spin wave solver
* [spinw.powspec](/spinw_powspec) powder spectrum calculator
* [spinw.energy](/spinw_energy) ground state energy
* [spinw.moment](/spinw_moment) moment reduction due to quantum fluctuations
* [spinw.spinwavesym](/spinw_spinwavesym) symbolic spin wave solver
* [spinw.symbolic](/spinw_symbolic) returns whether symbolic mode is on
* [spinw.meanfield](/spinw_meanfield) mean field calculation of q-dependent susceptibility (EXPERIMENTAL)
* [spinw.fourier](/spinw_fourier) fourier transformation of the Hamiltonian (EXPERIMENTAL)
* [spinw.fouriersym](/spinw_fouriersym) symbolic fourier transformation (EXPERIMENTAL)
 
#### Fitting spin wave spectrum

* [spinw.fitspec](/spinw_fitspec) fits spin wave energies using global optimization
* [spinw.matparser](/spinw_matparser) assigns new matrix values based on selectors
* [spinw.horace](/spinw_horace) outputs spectrum for Horace
    
#### Miscellaneous

* [spinw.copy](/spinw_copy) hard copy of SpinW object
* [spinw.export](/spinw_export) exports magnetic structure into file
* [spinw.fileid](/spinw_fileid) controls text output of commands
* [spinw.table](/spinw_table) formatted output of internal data
* [spinw.validate](/spinw_validate) validates SpinW object
* [spinw.version](/spinw_version) version of SpinW used to create the object
* [spinw.struct](/spinw_struct) convert SpinW ojbect to struct
* [spinw.clearcache](/spinw_clearcache) clear all data from cache, forcing recalculation (INTERNAL)
* [spinw.spinw](/spinw_spinw) constructor
 
 
#### Tutorials and documentation can be found here
#### <a href='https//www.psi.ch/spinw'>https//www.psi.ch/spinw</a>
#### Forum for questions
#### <a href='https//groups.google.com/forum/#!forum/spinwforum'>https//groups.google.com/forum/#!forum/spinwforum</a>
#### Lates version and bug reports/feature requests
#### <a href='https//github.com/tsdev/spinw'>https//github.com/tsdev/spinw</a>
 
  Reference page in Doc Center
     doc spinw

