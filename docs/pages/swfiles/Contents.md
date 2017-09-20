---
{title: Functions in swfiles, link: Functions in swfiles, summary: unstructured functions,
  keywords: sample, sidebar: sw_sidebar, permalink: swfiles.html, folder: swfiles,
  mathjax: 'true'}

---
This folder contains all the files related to SpinW but not yet split out
into separate libraries.
 
### Files
 
#### Transforming and plotting calculated spin wave spectrum
 
* [sw_econtract](sw_econtract.html) converts (Q,E) values to Q values for diffraction instrument
* [sw_egrid](sw_egrid.html) creates energy for spectrum color plot
* [sw_filelist](sw_filelist.html) lists spinw data in the Matlab workspace or in a .mat file
* [sw_instrument](sw_instrument.html) includes instrumental factors into the calculated spectrum
* [sw_magdomain](sw_magdomain.html) calculates the spin-spin correlation function for magnetic domains
* [sw_neutron](sw_neutron.html) calculates neutron scattering intensity for spin wave spectrum
* [sw_omegasum](sw_omegasum.html) removes degenerate and ghost magnon modes from spectrum
* [sw_plotspec](sw_plotspec.html) plots spin wave spectrum
* [sw_xray](sw_xray.html) calculates X-ray scattering intensity for phonon spectrum
 
#### Generate list of vectors in reciprocal space
 
* [sw_qgrid](sw_qgrid.html) creates a Q grid
* [sw_qscan](sw_qscan.html) creates continuous line between coordinates
 
#### Resolution claculation and convolution
 
* [sw_res](sw_res.html) reads a tabulated energy resolution from a file and fits with polynomial
* [sw_resconv](sw_resconv.html) convolution of a matrix and a Gaussian
* [sw_tofres](sw_tofres.html) includes Q resolution to the spectrum
 
#### SpinW model related functions
 
* [sw_extendlattice](sw_extendlattice.html) creates superlattice
* [sw_fstat](sw_fstat.html) calculates termodynamical averages during an annealing simulation
* [sw_intsf](sw_intsf.html) integrates the structure factor along given Q directions
* [sw_model](sw_model.html) creates different predefined spin models
* [sw_bonddim](sw_bonddim.html) find dimensionality of a periodic bond network
 
#### Constraint functions
 
* [gm_planar](gm_planar.html) planar magnetic structure constraint function 
* [gm_planard](gm_planard.html) planar magnetic structure constraint function 
* [gm_spherical3d](gm_spherical3d.html) magnetic structure constraint function with spherical parameterisation
* [gm_spherical3dd](gm_spherical3dd.html) magnetic structure constraint function with spherical parameterisation
 
#### Geometrical calculations
* [sw_angle](sw_angle.html) calculates the angle between 2 vectors
* [sw_cartesian](sw_cartesian.html) creates a right handed Cartesian coordinate system
* [sw_cmod](sw_cmod.html) modulo one with tolerance
* [sw_fsub](sw_fsub.html) simple graph vertex coloring
* [sw_mattype](sw_mattype.html) determines the type of square input matrix
* [sw_nvect](sw_nvect.html) determines the best normal vector for the set of vectors
* [sw_quadell](sw_quadell.html) calculates and plots the parameters of an ellipsoid from a quadratic form
* [sw_rot](sw_rot.html) rotates vectors around arbitrary axis in 3D
* [sw_rotmat](sw_rotmat.html) rotates vectors around arbitrary axis in 3D
* [sw_rotmatd](sw_rotmatd.html) rotates vectors around arbitrary axis in 3D
 
#### Text and graphical input/output for different high level commands
 
* [sw_annealfigure](sw_annealfigure.html) creates a figure for displaying the status of the annealing simulation
* [sw_annealplot](sw_annealplot.html) displays information about the annealing simulation
* [sw_label](sw_label.html) returns axis labels for spectrum plot
* [sw_circle](sw_circle.html) creates an array of the 3D coordinates of the circle circumference
* [sw_counter](sw_counter.html) print the number of calls to this function to the Command Line
* [sw_multicolor](sw_multicolor.html) creates RGB color data for multiple 2D overlapping plots
* [sw_parstr](sw_parstr.html) parses input string
* [sw_plotcell](sw_plotcell.html) plots cell structure with circles
* [sw_plotsf](sw_plotsf.html) plots the structure factor in the selected Q range in 1D or 2D
* [sw_status](sw_status.html) timer function that displays also the remaining time
 
#### Acessing the SpinW database
 
* [sw_atomdata](sw_atomdata.html) returns information on elements stored in the atom.dat file
* [sw_cff](sw_cff.html) returns the atomic charge form factor values for X-ray scattering
* [sw_mff](sw_mff.html) returns the magnetic form factor values and the coefficients
* [sw_nb](sw_nb.html) returns the bound coherent neutron scattering length (fm)
 
#### Symmetry calculations
 
* [sw_basismat](sw_basismat.html) determines allowed tensor components in a given point group symmetry
* [sw_mirror](sw_mirror.html) mirrors a 3D vector
 
#### Useful functions for physics
 
* [sw_bose](sw_bose.html) coefficient for boson correlation functions for different temperatures
* [sw_converter](sw_converter.html) converts energy and momentum units for a given particle
* [sw_fibo](sw_fibo.html) returns the last two Fibonacci number smaller or equal to the
 
#### Import functions
 
* [sw_import](sw_import.html) create SpinW object from .cif and FullProf Studio .fst files
* [sw_readspec](sw_readspec.html) read spin wave dispersion data from file
* [sw_readtable](sw_readtable.html) reads tabular data
 
#### Export functions
 
* [sw_idata](sw_idata.html) creates iData object
 
#### Miscellaneous
 
* [sw_freemem](sw_freemem.html) calculates the available memory
* [sw_initialize](sw_initialize.html) initializes spinw by removing user entries from the symmetry.dat file
* [sw_readparam](sw_readparam.html) parse input arguments (option, value pairs)
* [sw_rootdir](sw_rootdir.html) gives the path to the SpinW code
* [sw_uniquetol](sw_uniquetol.html) returns the unique column vectors within tolerance
* [sw_update](sw_update.html) updates the SpinW installation from the internet
* [sw_version](sw_version.html) returns the installed version of SpinW
* [sw_mex](sw_mex.html) compiles the mex files and test them
* [sw_notify](sw_notify.html) sends notification in OSX

