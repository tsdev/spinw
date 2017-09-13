---
{title: Class spinw, summary: SPINW class defines data structure and methods to calculate
    spin wave, keywords: sample, sidebar: sw_sidebar, permalink: spinw.html, folder: spinw,
  mathjax: 'true'}

---
SPINW class defines data structure and methods to calculate spin wave
dispersion in magnetic crystals.
 
obj = SPINW()
 
constructs a new [spinw](spinw.html) class object, with default parameters.
 
obj = SPINW(obj)
 
constructs new [spinw](spinw.html) class object. If obj is [spinw](spinw.html) class, it only
checks its data integrity. If obj is struct type, it creates new
[spinw](spinw.html) object and checks data integrity.
 
obj = SPINW(file_path)
 
construct new [spinw](spinw.html) class object, where file_path contains a string
of a .cif/.fst file defining an input crystal and/or magnetic
structure.
 
obj = SPINW(figure_handle)
 
copy the SpinW object stored in a previous structural plot.
 
The data structure behind the [spinw](spinw.html) object can be accessed by using
STRUCT(obj). All fields of the struct type data behind the [spinw](spinw.html)
object are accessible through the main field names of the obj object.
For example the lattice parameters:
  abc = obj.unit_cell.lat_const;
 
[spinw](spinw.html) is a handle class, that means that only the handle of the
object is copied in a swobj1 = swobj2 command. To create a copy
(clone) of an [spinw](spinw.html) object use:
   swobj1 = swobj2.COPY;
 
See also: SPINW.COPY, SPINW.STRUCT,
<a href='/Applications/MATLAB_R2012b.app/help/matlab/matlab_oop/comparing-handle-and-value-classes.html'>Comparing handle and value classes</a>.
 
### Methods
 
#### Lattice operations

* [[spinw](spinw.html).genlattice](/[spinw_genlattice](spinw_genlattice.html)) defines lattice parameters, angles and symmetry
* [[spinw](spinw.html).basisvector](/[spinw_basisvector](spinw_basisvector.html)) returns the basis vectors in a matrix
* [[spinw](spinw.html).rl](/[spinw_rl](spinw_rl.html)) returns the reciprocal lattice vectors in a matrix
* [[spinw](spinw.html).nosym](/[spinw_nosym](spinw_nosym.html)) reduces the lattice symmetry to P0
* [[spinw](spinw.html).newcell](/[spinw_newcell](spinw_newcell.html)) generates an arbitrary new cell and reduces symmetry to P0
* [[spinw](spinw.html).addatom](/[spinw_addatom](spinw_addatom.html)) adds a new atoms to the lattice
* [[spinw](spinw.html).unitcell](/[spinw_unitcell](spinw_unitcell.html)) selects a subset of the symmetry inequivalent atoms
* [[spinw](spinw.html).abc](/[spinw_abc](spinw_abc.html)) returns a vector [a b c alpha beta gamma]
* [[spinw](spinw.html).atom](/[spinw_atom](spinw_atom.html)) returns all symmetry generated atom
* [[spinw](spinw.html).matom](/[spinw_matom](spinw_matom.html)) returns all symmetry generated magnetic atoms
* [[spinw](spinw.html).natom](/[spinw_natom](spinw_natom.html)) returns the number of symmetry inequivalent atoms
* [[spinw](spinw.html).formula](/[spinw_formula](spinw_formula.html)) returns basic unit cell information
* [[spinw](spinw.html).disp](/[spinw_disp](spinw_disp.html)) return basic lattice information
* [[spinw](spinw.html).symmetry](/[spinw_symmetry](spinw_symmetry.html)) returns whether symmetry is on (>P0)
    
#### Plotting

* [[spinw](spinw.html).plot](/[spinw_plot](spinw_plot.html)) plots crystal, magnetic structure and spin Hamiltonian
 
#### Crystallographic twin operations

* [[spinw](spinw.html).addtwin](/[spinw_addtwin](spinw_addtwin.html)) adds crystallographic twin
* [[spinw](spinw.html).twinq](/[spinw_twinq](spinw_twinq.html)) determines reciprocal lattice positions in twins
* [[spinw](spinw.html).notwin](/[spinw_notwin](spinw_notwin.html)) removes the twins
* [[spinw](spinw.html).ntwin](/[spinw_ntwin](spinw_ntwin.html)) number of twins
 
#### Magnetic structure operations

* [[spinw](spinw.html).genmagstr](/[spinw_genmagstr](spinw_genmagstr.html)) generate different types of magnetic structures
* [[spinw](spinw.html).magstr](/[spinw_magstr](spinw_magstr.html)) returns the magnetic structure in a rotating frame representation
* [[spinw](spinw.html).magtable](/[spinw_magtable](spinw_magtable.html)) returns a table of the magnetic structure
* [[spinw](spinw.html).nmagext](/[spinw_nmagext](spinw_nmagext.html)) returns the number of magnetic moments in the supercell
* [[spinw](spinw.html).optmagstr](/[spinw_optmagstr](spinw_optmagstr.html)) optimizes magnetic structure with constraints
* [[spinw](spinw.html).optmagk](/[spinw_optmagk](spinw_optmagk.html)) optimizes the magnetic propagation vector
* [[spinw](spinw.html).optmagsteep](/[spinw_optmagsteep](spinw_optmagsteep.html)) optimizes the moment directions for fixed-k
* [[spinw](spinw.html).anneal](/[spinw_anneal](spinw_anneal.html)) simulated annealing of magnetic structures
* [[spinw](spinw.html).annealloop](/[spinw_annealloop](spinw_annealloop.html)) simulated annealing and scanning a parameter in the Hamiltonian
* [[spinw](spinw.html).structfact](/[spinw_structfact](spinw_structfact.html)) calculates magnetic structure factor (EXPERIMENTAL)
    
#### Matrix operations

* [[spinw](spinw.html).addmatrix](/[spinw_addmatrix](spinw_addmatrix.html)) add a new matrix to the internal data
* [[spinw](spinw.html).getmatrix](/[spinw_getmatrix](spinw_getmatrix.html)) determines the symmetry allowed elements of a matrix
* [[spinw](spinw.html).setmatrix](/[spinw_setmatrix](spinw_setmatrix.html)) generates matrix values based on symmetry allowed elements
* [[spinw](spinw.html).nmat](/[spinw_nmat](spinw_nmat.html)) number of matrices
    
#### Spin Hamiltonian generations

* [[spinw](spinw.html).quickham](/[spinw_quickham](spinw_quickham.html)) quick generation of Heisenberg Hamiltonian
* [[spinw](spinw.html).gencoupling](/[spinw_gencoupling](spinw_gencoupling.html)) generates the list of bonds
* [[spinw](spinw.html).addcoupling](/[spinw_addcoupling](spinw_addcoupling.html)) assigns a matrix to a bond
* [[spinw](spinw.html).couplingtable](/[spinw_couplingtable](spinw_couplingtable.html)) lists information on the bonds
* [[spinw](spinw.html).addaniso](/[spinw_addaniso](spinw_addaniso.html)) assigns a matrix to a magnetic ion as anisotropy
* [[spinw](spinw.html).addg](/[spinw_addg](spinw_addg.html)) assigns a matrix to a magnetic ion as g-tensor
* [[spinw](spinw.html).field](/[spinw_field](spinw_field.html)) stores the magnetic field
* [[spinw](spinw.html).nbond](/[spinw_nbond](spinw_nbond.html)) number of bonds
* [[spinw](spinw.html).temperature](/[spinw_temperature](spinw_temperature.html)) temperature for thermal population calculation
* [[spinw](spinw.html).intmatrix](/[spinw_intmatrix](spinw_intmatrix.html)) returns the spin Hamiltonian after symmetry applied (INTERNAL)
* [[spinw](spinw.html).symop](/[spinw_symop](spinw_symop.html)) generates the symmetry operators on bonds and magnetic atoms (INTERNAL)
* [[spinw](spinw.html).setunit](/[spinw_setunit](spinw_setunit.html)) sets the physical units
    
#### Calculators

* [[spinw](spinw.html).spinwave](/[spinw_spinwave](spinw_spinwave.html)) linear spin wave solver
* [[spinw](spinw.html).powspec](/[spinw_powspec](spinw_powspec.html)) powder spectrum calculator
* [[spinw](spinw.html).energy](/[spinw_energy](spinw_energy.html)) ground state energy
* [[spinw](spinw.html).moment](/[spinw_moment](spinw_moment.html)) moment reduction due to quantum fluctuations
* [[spinw](spinw.html).spinwavesym](/[spinw_spinwavesym](spinw_spinwavesym.html)) symbolic spin wave solver
* [[spinw](spinw.html).symbolic](/[spinw_symbolic](spinw_symbolic.html)) returns whether symbolic mode is on
* [[spinw](spinw.html).meanfield](/[spinw_meanfield](spinw_meanfield.html)) mean field calculation of q-dependent susceptibility (EXPERIMENTAL)
* [[spinw](spinw.html).fourier](/[spinw_fourier](spinw_fourier.html)) fourier transformation of the Hamiltonian (EXPERIMENTAL)
* [[spinw](spinw.html).fouriersym](/[spinw_fouriersym](spinw_fouriersym.html)) symbolic fourier transformation (EXPERIMENTAL)
 
#### Fitting spin wave spectrum

* [[spinw](spinw.html).fitspec](/[spinw_fitspec](spinw_fitspec.html)) fits spin wave energies using global optimization
* [[spinw](spinw.html).matparser](/[spinw_matparser](spinw_matparser.html)) assigns new matrix values based on selectors
* [[spinw](spinw.html).horace](/[spinw_horace](spinw_horace.html)) outputs spectrum for Horace
    
#### Miscellaneous

* [[spinw](spinw.html).copy](/[spinw_copy](spinw_copy.html)) hard copy of SpinW object
* [[spinw](spinw.html).export](/[spinw_export](spinw_export.html)) exports magnetic structure into file
* [[spinw](spinw.html).fileid](/[spinw_fileid](spinw_fileid.html)) controls text output of commands
* [[spinw](spinw.html).table](/[spinw_table](spinw_table.html)) formatted output of internal data
* [[spinw](spinw.html).validate](/[spinw_validate](spinw_validate.html)) validates SpinW object
* [[spinw](spinw.html).version](/[spinw_version](spinw_version.html)) version of SpinW used to create the object
* [[spinw](spinw.html).struct](/[spinw_struct](spinw_struct.html)) convert SpinW ojbect to struct
* [[spinw](spinw.html).clearcache](/spinw_clearcache) clear all data from cache, forcing recalculation (INTERNAL)
* [[spinw](spinw.html).spinw](/[spinw_spinw](spinw_spinw.html)) constructor
 
 
#### Tutorials and documentation can be found here
#### <a href='https//www.psi.ch/[spinw](spinw.html)'>https//www.psi.ch/[spinw](spinw.html)</a>
#### Forum for questions
#### <a href='https//groups.google.com/forum/#!forum/spinwforum'>https//groups.google.com/forum/#!forum/spinwforum</a>
#### Lates version and bug reports/feature requests
#### <a href='https//github.com/tsdev/[spinw](spinw.html)'>https//github.com/tsdev/[spinw](spinw.html)</a>
 
  Reference page in Doc Center
     doc [spinw](spinw.html)

