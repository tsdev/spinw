---
{title: spinw( ), keywords: sample, summary: SPINW class defines data structure and methods to calculate spin wave,
  sidebar: sw_sidebar, permalink: '@spinw_spinw.html', folder: '@spinw', mathjax: 'true'}

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
 
  spinw Methods:
  Lattice operations:
    genlattice      - defines lattice parameters, angles and symmetry
    basisvector     - returns the basis vectors in a matrix
    rl              - returns the reciprocal lattice vectors in a matrix
    nosym           - reduces the lattice symmetry to P0
    newcell         - generates an arbitrary new cell and reduces symmetry to P0
    addatom         - adds a new atoms to the lattice
    unitcell        - selects a subset of the symmetry inequivalent atoms 
    abc             - returns a vector [a b c alpha beta gamma]
    atom            - returns all symmetry generated atom
    matom           - returns all symmetry generated magnetic atoms
    natom           - returns the number of symmetry inequivalent atoms
    formula         - returns basic unit cell information
    disp            - return basic lattice information
    symmetry        - returns whether symmetry is on (>P0)
    
  Plotting:
    plot            - plots crystal, magnetic structure and spin Hamiltonian
 
  Crystallographic twin operations:
    addtwin         - adds crystallographic twin
    twinq           - determines reciprocal lattice positions in twins
    notwin          - removes the twins
    ntwin           - number of twins
 
  Magnetic structure operations:
    genmagstr       - generate different types of magnetic structures
    magstr          - returns the magnetic structure in a rotating frame representation
    magtable        - returns a table of the magnetic structure
    nmagext         - returns the number of magnetic moments in the supercell
    optmagstr       - optimizes magnetic structure with constraints
    optmagk         - optimizes the magnetic propagation vector
    optmagsteep     - optimizes the moment directions for fixed-k
    anneal          - simulated annealing of magnetic structures
    annealloop      - simulated annealing and scanning a parameter in the Hamiltonian
    structfact      - calculates magnetic structure factor (EXPERIMENTAL)
    
  Matrix operations:
    addmatrix        - add a new matrix to the internal data
    getmatrix       - determines the symmetry allowed elements of a matrix
    setmatrix       - generates matrix values based on symmetry allowed elements
    nmat            - number of matrices
    
  Spin Hamiltonian generations:
    quickham        - quick generation of Heisenberg Hamiltonian
    gencoupling     - generates the list of bonds
    addcoupling     - assigns a matrix to a bond
    couplingtable   - lists information on the bonds
    addaniso        - assigns a matrix to a magnetic ion as anisotropy
    addg            - assigns a matrix to a magnetic ion as g-tensor
    field           - stores the magnetic field
    nbond           - number of bonds
    temperature     - temperature for thermal population calculation
    intmatrix       - returns the spin Hamiltonian after symmetry applied (INTERNAL)
    symop           - generates the symmetry operators on bonds and magnetic atoms (INTERNAL)
    
  Calculators:
    spinwave        - linear spin wave solver
    powspec         - powder spectrum calculator
    energy          - ground state energy
    moment          - moment reduction due to quantum fluctuations
    spinwavesym     - symbolic spin wave solver
    symbolic        - returns whether symbolic mode is on
    meanfield       - mean field calculation of q-dependent susceptibility (EXPERIMENTAL)
    fourier         - fourier transformation of the Hamiltonian (EXPERIMENTAL)
    fouriersym      - symbolic fourier transformation (EXPERIMENTAL)
 
  Fitting spin wave spectrum:
    fitspec         - fits spin wave energies using global optimization
    matparser       - assigns new matrix values based on selectors
    horace          - outputs spectrum for Horace
    
  Miscellaneous:
    copy            - hard copy of SpinW object
    export          - exports magnetic structure into file
    fileid          - controls text output of commands
    table           - formatted output of internal data
    validate        - validates SpinW object
    version         - version of SpinW used to create the object
    struct          - convert SpinW ojbect to struct
    clearcache      - clear all data from cache, forcing recalculation (INTERNAL)
    spinw           - constructor
 
 
  Tutorials and documentation can be found here:
  <a href='https://www.psi.ch/spinw'>https://www.psi.ch/spinw</a>
  Forum for questions:
  <a href='https://groups.google.com/forum/#!forum/spinwforum'>https://groups.google.com/forum/#!forum/spinwforum</a>
  Lates version and bug reports/feature requests:
  <a href='https://github.com/tsdev/spinw'>https://github.com/tsdev/spinw</a>
 

    Reference page in Doc Center
       doc spinw

