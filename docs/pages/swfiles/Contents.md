---
{title: Functions in swfiles, link: Functions in swfiles, summary: unstructured functions,
  keywords: sample, sidebar: sw_sidebar, permalink: swfiles, folder: swfiles, mathjax: 'true'}

---
This folder contains all the files related to SpinW but not yet split out
into separate libraries.
 
### Files
 
#### Transforming and plotting calculated spin wave spectrum
 
* [sw_econtract](sw_econtract) converts (Q,E) values to Q values for diffraction instrument
* [sw_egrid](sw_egrid) calculates energy bins of a spectrum 
* [sw_filelist](sw_filelist) lists spinw data in the Matlab workspace or in a .mat file
* [sw_instrument](sw_instrument) convolutes spectrum with different functions
* [sw_magdomain](sw_magdomain) calculates the spin-spin correlation function for magnetic domains
* [sw_neutron](sw_neutron) calculates neutron scattering intensity for spin wave spectrum
* [sw_omegasum](sw_omegasum) removes degenerate and ghost magnon modes from spectrum
* [sw_plotspec](sw_plotspec) plots spectrum
* [sw_xray](sw_xray) calculates X-ray scattering intensity for phonon spectrum
 
#### Generate list of vectors in reciprocal space
 
* [sw_qgrid](sw_qgrid) creates a Q grid
* [sw_qscan](sw_qscan) creates continuous line between coordinates
 
#### Resolution claculation and convolution
 
* [sw_res](sw_res) reads a tabulated energy resolution from a file and fits with polynomial
* [sw_resconv](sw_resconv) convolution of a matrix and a Gaussian
* [sw_tofres](sw_tofres) includes Q resolution to the spectrum
 
#### SpinW model related functions
 
* [sw_extendlattice](sw_extendlattice) creates superlattice
* [sw_fstat](sw_fstat) calculates termodynamical averages during an annealing simulation
* [sw_intsf](sw_intsf) integrates the structure factor along given Q directions
* [sw_model](sw_model) creates different predefined spin models
* [sw_bonddim](sw_bonddim) find dimensionality of a periodic bond network
 
#### Constraint functions
 
* [gm_planar](gm_planar) planar magnetic structure constraint function 
* [gm_planard](gm_planard) planar magnetic structure constraint function 
* [gm_spherical3d](gm_spherical3d) magnetic structure constraint function with spherical parameterisation
* [gm_spherical3dd](gm_spherical3dd) magnetic structure constraint function with spherical parameterisation
 
#### Geometrical calculations
* [sw_angle](sw_angle) calculates the angle between 2 vectors
* [sw_cartesian](sw_cartesian) creates a right handed Cartesian coordinate system
* [sw_cmod](sw_cmod) modulo one with tolerance
* [sw_fsub](sw_fsub) simple graph vertex coloring
* [sw_mattype](sw_mattype) determines the type of square input matrix
* [sw_nvect](sw_nvect) determines the best normal vector for the set of vectors
* [sw_quadell](sw_quadell) calculates and plots the parameters of an ellipsoid from a quadratic form
* [sw_rot](sw_rot) rotates vectors around arbitrary axis in 3D
* [sw_rotmat](sw_rotmat) rotates vectors around arbitrary axis in 3D
* [sw_rotmatd](sw_rotmatd) rotates vectors around arbitrary axis in 3D
 
#### Text and graphical input/output for different high level commands
 
* [sw_annealfigure](sw_annealfigure) creates a figure for displaying the status of the annealing simulation
* [sw_annealplot](sw_annealplot) displays information about the annealing simulation
* [sw_label](sw_label) returns axis labels for spectrum plot
* [sw_circle](sw_circle) creates an array of the 3D coordinates of the circle circumference
* [sw_counter](sw_counter) print the number of calls to this function to the Command Line
* [sw_multicolor](sw_multicolor) creates RGB color data for multiple 2D overlapping plots
* [sw_parstr](sw_parstr) parses input string
* [sw_plotcell](sw_plotcell) plots cell structure with circles
* [sw_plotsf](sw_plotsf) plots the structure factor in the selected Q range in 1D or 2D
* [sw_status](sw_status) timer and remaining time estimator
 
#### Acessing the SpinW database
 
* [sw_atomdata](sw_atomdata) returns information on elements stored in the atom.dat file
* [sw_cff](sw_cff) returns the atomic charge form factor values for X-ray scattering
* [sw_mff](sw_mff) returns the magnetic form factor values and the coefficients
* [sw_nb](sw_nb) returns the bound coherent neutron scattering length (fm)
 
#### Symmetry calculations
 
* [sw_basismat](sw_basismat) determines allowed tensor components in a given point group symmetry
* [sw_mirror](sw_mirror) mirrors a 3D vector
 
#### Useful functions for physics
 
* [sw_bose](sw_bose) coefficient for boson correlation functions for different temperatures
* [sw_converter](sw_converter) converts energy and momentum units for a given particle
* [sw_fibo](sw_fibo) returns the last two Fibonacci number smaller or equal to the
 
#### Import functions
 
* [sw_import](sw_import) create SpinW object from .cif and FullProf Studio .fst files
* [sw_readspec](sw_readspec) read spin wave dispersion data from file
* [sw_readtable](sw_readtable) reads tabular data
 
#### Export functions
 
* [sw_idata](sw_idata) creates iData object
 
#### Miscellaneous
 
* [sw_freemem](sw_freemem) calculates the available memory
* [sw_initialize](sw_initialize) initializes spinw by removing user entries from the symmetry.dat file
* [sw_readparam](sw_readparam) parse input arguments (option, value pairs)
* [sw_rootdir](sw_rootdir) gives the path to the SpinW code
* [sw_uniquetol](sw_uniquetol) returns the unique column vectors within tolerance
* [sw_update](sw_update) updates the SpinW installation from the internet
* [sw_version](sw_version) returns the installed version of SpinW
* [sw_mex](sw_mex) compiles the mex files and test them
* [sw_notify](sw_notify) sends notification in OSX

{% include links.html %}
