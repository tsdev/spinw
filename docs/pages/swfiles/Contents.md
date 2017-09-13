---
{title: Functions in swfiles, summary: This folder contains all the files related
    to SpinW but not yet split out, keywords: sample, sidebar: sw_sidebar, permalink: swfiles.html,
  folder: swfiles, mathjax: 'true'}

---
This folder contains all the files related to SpinW but not yet split out
into separate libraries.
 
### Files
 
#### Transforming and plotting calculated spin wave spectrum
 

* [[sw_econtract](sw_econtract.html)](/sw_econtract) converts (Q,omega) values to Qm values for diffraction instrument
* [[sw_egrid](sw_egrid.html)](/sw_egrid) creates energy for spectrum color plot
* [[sw_filelist](sw_filelist.html)](/sw_filelist) lists [spinw](spinw.html) data in the Matlab workspace or in a .mat file
* [[sw_instrument](sw_instrument.html)](/sw_instrument) includes instrumental factors into the calculated spectrum
* [[sw_magdomain](sw_magdomain.html)](/sw_magdomain) calculates the spin-spin correlation function for magnetic domains
* [[sw_neutron](sw_neutron.html)](/sw_neutron) calculates neutron scattering intensity for spin wave spectrum
* [[sw_omegasum](sw_omegasum.html)](/sw_omegasum) removes degenerate and ghost magnon modes from spectrum
* [[sw_plotspec](sw_plotspec.html)](/sw_plotspec) plots spin wave spectrum
* [[sw_xray](sw_xray.html)](/sw_xray) calculates X-ray scattering intensity for phonon spectrum
 
#### Generate list of vectors in reciprocal space
 

* [[sw_qgrid](sw_qgrid.html)](/sw_qgrid) creates a Q grid
* [[sw_qscan](sw_qscan.html)](/sw_qscan) creates linear scans between Q points in 3D
 
#### Resolution claculation and convolution
 

* [[sw_res](sw_res.html)](/sw_res) reads a tabulated energy resolution from a file and fits with polynomial
* [[sw_resconv](sw_resconv.html)](/sw_resconv) Convolute Gaussian with variable width along the first dimension of a matrix
* [[sw_tofres](sw_tofres.html)](/sw_tofres) includes Q resolution to the spectrum
 
#### SpinW model related functions
 

* [[sw_extendlattice](sw_extendlattice.html)](/sw_extendlattice) creates superlattice
* [[sw_fstat](sw_fstat.html)](/sw_fstat) calculates termodynamical averages during an annealing simulation
* [[sw_intsf](sw_intsf.html)](/sw_intsf) integrates the structure factor along given Q directions
* [[sw_model](sw_model.html)](/sw_model) creates different predefined spin models
* [[sw_bonddim](sw_bonddim.html)](/sw_bonddim) find dimensionality of a periodic bond network
 
#### Constraint functions
 

* [[gm_planar](gm_planar.html)](/gm_planar) planar magnetic structure constraint function
* [[gm_planard](gm_planard.html)](/gm_planard) planar magnetic structure constraint function
* [[gm_spherical3d](gm_spherical3d.html)](/gm_spherical3d) magnetic structure constraint function with spherical parameterisation
* [[gm_spherical3dd](gm_spherical3dd.html)](/gm_spherical3dd) magnetic structure constraint function with spherical parameterisation
 
#### Geometrical calculations

* [[sw_angle](sw_angle.html)](/sw_angle) calculates the angle between 2 vectors
* [[sw_cartesian](sw_cartesian.html)](/sw_cartesian) creates a right handed Cartesian coordinate system
* [[sw_cmod](sw_cmod.html)](/sw_cmod) modulo one with tolerance
* [[sw_fsub](sw_fsub.html)](/sw_fsub) simple graph vertex coloring
* [[sw_mattype](sw_mattype.html)](/sw_mattype) determines the type of square input matrix
* [[sw_nvect](sw_nvect.html)](/sw_nvect) determines the best normal vector for the set of vectors
* [[sw_quadell](sw_quadell.html)](/sw_quadell) calculates and plots the parameters of an ellipsoid from a quadratic form
* [[sw_rot](sw_rot.html)](/sw_rot) rotates vectors around arbitrary axis in 3D
* [[sw_rotmat](sw_rotmat.html)](/sw_rotmat) rotates vectors around arbitrary axis in 3D
* [[sw_rotmatd](sw_rotmatd.html)](/sw_rotmatd) rotates vectors around arbitrary axis in 3D
 
#### Text and graphical input/output for different high level commands
 

* [[sw_annealfigure](sw_annealfigure.html)](/sw_annealfigure) creates a figure for displaying the status of the annealing simulation
* [[sw_annealplot](sw_annealplot.html)](/sw_annealplot) displays information about the annealing simulation
* [[sw_label](sw_label.html)](/sw_label) returns axis labels for spectrum plot
* [[sw_circle](sw_circle.html)](/sw_circle) creates an array of the 3D coordinates of the circle circumference
* [[sw_counter](sw_counter.html)](/sw_counter) print the number of calls to this functio to the Command Line
* [[sw_multicolor](sw_multicolor.html)](/sw_multicolor) creates RGB color data for multiple 2D overlapping plots
* [[sw_parstr](sw_parstr.html)](/sw_parstr) parses input string
* [[sw_plotcell](sw_plotcell.html)](/sw_plotcell) plots cell structure with circles
* [[sw_plotsf](sw_plotsf.html)](/sw_plotsf) plots the structure factor in the selected Q range in 1D or 2D
* [[sw_status](sw_status.html)](/sw_status) timer function that displays also the remaining time
 
#### Acessing the SpinW database
 

* [[sw_atomdata](sw_atomdata.html)](/sw_atomdata) returns information on elements stored in the atom.dat file
* [[sw_cff](sw_cff.html)](/sw_cff) returns the atomic charge form factor values for X-ray scattering
* [[sw_mff](sw_mff.html)](/sw_mff) returns the magnetic form factor values and the coefficients
* [[sw_nb](sw_nb.html)](/sw_nb) returns the bound coherent neutron scattering length (fm)
 
#### Symmetry calculations
 

* [[sw_basismat](sw_basismat.html)](/sw_basismat) determines allowed tensor components in a given point group symmetry
* [[sw_mirror](sw_mirror.html)](/sw_mirror) mirrors a 3D vector
 
#### Useful functions for physics
 

* [[sw_bose](sw_bose.html)](/sw_bose) coefficient for boson correlation functions for different temperatures
* [[sw_converter](sw_converter.html)](/sw_converter) converts energy and momentum units for a given particle
* [[sw_fibo](sw_fibo.html)](/sw_fibo) returns the last two Fibonacci number smaller or equal to the
 
#### Import functions
 

* [[sw_import](sw_import.html)](/sw_import) create SpinW object from .cif and FullProf Studio .fst files
* [[sw_readspec](sw_readspec.html)](/sw_readspec) read spin wave dispersion data from file
* [[sw_readtable](sw_readtable.html)](/sw_readtable) reads tabular data
 
#### Export functions
 

* [[sw_idata](sw_idata.html)](/sw_idata) creates iData object
 
#### Miscellaneous
 

* [[sw_freemem](sw_freemem.html)](/sw_freemem) gives the amount of free RAM in bytes
* [[sw_initialize](sw_initialize.html)](/sw_initialize) initializes [spinw](spinw.html) by removing user entries from the symmetry.dat file
* [[sw_readparam](sw_readparam.html)](/sw_readparam) parse input arguments (option, value pairs)
* [[sw_rootdir](sw_rootdir.html)](/sw_rootdir) gives the path to the SpinW code
* [[sw_uniquetol](sw_uniquetol.html)](/sw_uniquetol) returns the unique column vectors within tolerance
* [[sw_update](sw_update.html)](/sw_update) updates the SpinW installation from the internet
* [[sw_version](sw_version.html)](/sw_version) returns the installed version of SpinW
* [[sw_mex](sw_mex.html)](/sw_mex) compiles the mex files and test them
* [[sw_notify](sw_notify.html)](/sw_notify) sends system notification in OSX

