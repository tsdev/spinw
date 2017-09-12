---
{title: Functions in swfiles, summary: This folder contains all the files related
    to SpinW but not yet split out, keywords: sample, sidebar: sw_sidebar, permalink: swfiles.html,
  folder: swfiles, mathjax: 'true'}

---
This folder contains all the files related to SpinW but not yet split out
into separate libraries.
 
### Files
 
#### Transforming and plotting calculated spin wave spectrum:
 

* [sw_econtract](/swfiles_sw_econtract) converts (Q,omega) values to Qm values for diffraction instrument
* [sw_egrid](/swfiles_sw_egrid) creates energy for spectrum color plot
* [sw_filelist](/swfiles_sw_filelist) lists spinw data in the Matlab workspace or in a .mat file
* [sw_instrument](/swfiles_sw_instrument) includes instrumental factors into the calculated spectrum
  sw_magdomain     - calculates the spin-spin correlation function for magnetic domains
* [sw_neutron](/swfiles_sw_neutron) calculates neutron scattering intensity for spin wave spectrum
* [sw_omegasum](/swfiles_sw_omegasum) removes degenerate and ghost magnon modes from spectrum
* [sw_plotspec](/swfiles_sw_plotspec) plots spin wave spectrum
  sw_xray          - calculates X-ray scattering intensity for phonon spectrum
 
#### Generate list of vectors in reciprocal space:
 

* [sw_qgrid](/swfiles_sw_qgrid) creates a Q grid
* [sw_qscan](/swfiles_sw_qscan) creates linear scans between Q points in 3D
 
#### Resolution claculation and convolution:
 

* [sw_res](/swfiles_sw_res) reads a tabulated energy resolution from a file and fits with polynomial
* [sw_resconv](/swfiles_sw_resconv) Convolute Gaussian with variable width along the first dimension of a matrix
* [sw_tofres](/swfiles_sw_tofres) includes Q resolution to the spectrum
 
#### SpinW model related functions:
 

* [sw_extendlattice](/swfiles_sw_extendlattice) creates superlattice
* [sw_fstat](/swfiles_sw_fstat) calculates termodynamical averages during an annealing simulation
* [sw_intsf](/swfiles_sw_intsf) integrates the structure factor along given Q directions
* [sw_model](/swfiles_sw_model) creates different predefined spin models
* [sw_bonddim](/swfiles_sw_bonddim) find dimensionality of a periodic bond network
 
#### Constraint functions for spinw.optmagstr():
 

* [gm_planar](/swfiles_gm_planar) planar magnetic structure constraint function
* [gm_planard](/swfiles_gm_planard) planar magnetic structure constraint function
* [gm_spherical3d](/swfiles_gm_spherical3d) magnetic structure constraint function with spherical parameterisation
* [gm_spherical3dd](/swfiles_gm_spherical3dd) magnetic structure constraint function with spherical parameterisation
 
#### Geometrical calculations:

* [sw_angle](/swfiles_sw_angle) calculates the angle between 2 vectors
* [sw_cartesian](/swfiles_sw_cartesian) creates a right handed Cartesian coordinate system
* [sw_cmod](/swfiles_sw_cmod) modulo one with tolerance
* [sw_fsub](/swfiles_sw_fsub) simple graph vertex coloring
* [sw_mattype](/swfiles_sw_mattype) determines the type of square input matrix
* [sw_nvect](/swfiles_sw_nvect) determines the best normal vector for the set of vectors
* [sw_quadell](/swfiles_sw_quadell) calculates and plots the parameters of an ellipsoid from a quadratic form
* [sw_rot](/swfiles_sw_rot) rotates vectors around arbitrary axis in 3D
* [sw_rotmat](/swfiles_sw_rotmat) rotates vectors around arbitrary axis in 3D
* [sw_rotmatd](/swfiles_sw_rotmatd) rotates vectors around arbitrary axis in 3D
 
#### Text and graphical input/output for different high level commands:
 

* [sw_annealfigure](/swfiles_sw_annealfigure) creates a figure for displaying the status of the annealing simulation
* [sw_annealplot](/swfiles_sw_annealplot) displays information about the annealing simulation
* [sw_label](/swfiles_sw_label) returns axis labels for spectrum plot
* [sw_circle](/swfiles_sw_circle) creates an array of the 3D coordinates of the circle circumference
* [sw_counter](/swfiles_sw_counter) print the number of calls to this functio to the Command Line
* [sw_multicolor](/swfiles_sw_multicolor) creates RGB color data for multiple 2D overlapping plots
* [sw_parstr](/swfiles_sw_parstr) parses input string
* [sw_plotcell](/swfiles_sw_plotcell) plots cell structure with circles
* [sw_plotsf](/swfiles_sw_plotsf) plots the structure factor in the selected Q range in 1D or 2D
* [sw_status](/swfiles_sw_status) timer function that displays also the remaining time
 
#### Acessing the SpinW database:
 

* [sw_atomdata](/swfiles_sw_atomdata) returns information on elements stored in the atom.dat file
  sw_cff           - returns the atomic charge form factor values for X-ray scattering
* [sw_mff](/swfiles_sw_mff) returns the magnetic form factor values and the coefficients
* [sw_nb](/swfiles_sw_nb) returns the bound coherent neutron scattering length (fm)
 
#### Symmetry calculations:
 

* [sw_basismat](/swfiles_sw_basismat) determines allowed tensor components in a given point group symmetry
* [sw_mirror](/swfiles_sw_mirror) mirrors a 3D vector
 
#### Useful functions for physics:
 

* [sw_bose](/swfiles_sw_bose) coefficient for boson correlation functions for different temperatures
* [sw_converter](/swfiles_sw_converter) converts energy and momentum units for a given particle
* [sw_fibo](/swfiles_sw_fibo) returns the last two Fibonacci number smaller or equal to the
 
#### Import functions:
 

* [sw_import](/swfiles_sw_import) create SpinW object from .cif and FullProf Studio .fst files
* [sw_readspec](/swfiles_sw_readspec) read spin wave dispersion data from file
* [sw_readtable](/swfiles_sw_readtable) reads tabular data
 
#### Export functions:
 

* [sw_idata](/swfiles_sw_idata) creates iData object
 
#### Other files:
 

* [sw_freemem](/swfiles_sw_freemem) gives the amount of free RAM in bytes
* [sw_initialize](/swfiles_sw_initialize) initializes spinw by removing user entries from the symmetry.dat file
* [sw_readparam](/swfiles_sw_readparam) parse input arguments (option, value pairs)
* [sw_rootdir](/swfiles_sw_rootdir) gives the path to the SpinW code
* [sw_uniquetol](/swfiles_sw_uniquetol) returns the unique column vectors within tolerance
* [sw_update](/swfiles_sw_update) updates the SpinW installation from the internet
* [sw_version](/swfiles_sw_version) returns the installed version of SpinW
* [sw_mex](/swfiles_sw_mex) compiles the mex files and test them
* [sw_notify](/swfiles_sw_notify) sends system notification in OSX
