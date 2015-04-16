Introduction to SpinW
=====================

SpinW is a Matlab library that can solve spin Hamiltonians using mean field theory and linear spin wave theory. It is optimised for spin wave calculation on complex lattices that one often encounter experimentally.


Model
--------
In short SpinW can solve the following spin Hamiltonian using classical and quasi classical numerical methods: 

..  math::
    \mathcal{H}=\sum_{i,j}\mathbf{s}^\intercal_i\cdot J_{ij}\cdot \mathbf{s}_j + \sum_i \mathbf{s}^\intercal_i\cdot A_i\cdot \mathbf{s}_i + \mu_B\mathbf{B}^\intercal\cdot \sum_i\mathbf{g}_i\cdot \mathbf{s}_i

where :math:`\mathbf{s}_i` are spin vector operators, :math:`J_{ij}` are 3x3 exchange matrices describing coupling between spins, :math:`A_{ij}` are 3x3 single ion anisotropy matrices, :math:`\mathbf{B}` is external magnetic field, :math:`g_i` are the g-tensors with 3x3 elements and :math:`\mu_B` is the Bohr-magneton. The 3x3 matrices are necessary to include all possible anisotropic and antisymmetric interactions.


Features
--------

Crystal structures
....................

* definition of crystal lattice with arbitrary unit cell, using symmetry operators
* definition of non-magnetic atoms and magnetic atoms with arbitrary moment size
* publication quality plotting of crystal structures (atoms, labels, axes, surrounding polyhedron, anisotropy ellipsoids, DM vector, etc.)

Magnetic structures
...................
* definition of 1D, 2D and 3D magnetic structures
* representation of single-Q incommensurate structures using rotating coordinate system or complex basis vectors
* generation of magnetic structures on a magnetic supercell
* plotting of magnetic structures
Magnetic interactions
......................
* simple assignment of magnetic interactions to bonds based on bond length
* possible interactions: Heisenberg, Dzyaloshinskii-Moriya, anisotropic and general 3x3 exchange tensor
* arbitrary single ion anisotropy tensor (easy-plane, easy-axis, etc.)
* Zeeman energy in homogeneous magnetic field including arbitrary g-tensor
* calculation of symmetry allowed elements of the above tensors based on the crystallographic space group
Magnetic ground state optimisation
..................................
* minimization of the classical energy  assuming single-Q magnetic structure
* simulated annealing using the Metropolis algorithm on large magnetic supercells
* calculation of properties in thermodynamical equilibrium (heat capacity, magnetic susceptibility, etc.)
* magnetic structure factor calculation using FFT
* simulation of magnetic neutron diffraction and diffuse scattering
Spin wave simulation
..............................................................................
* solution for sommensurate and single-Q incommensurate magnetic structures
* calculation of spin wave dispersion, spin-spin correlation functions
* calculation of neutron scattering cross section for unpolarized neutrons including the magnetic form factor
* calculation of polarized neutron scattering cross sections
* possibility to include different moment sizes for different magnetic atoms
* calculation of powder averaged spin wave spectrum
Plotting of spin wave spectrum
...............................
* plotting of dispersions and correlation functions
* calculation and plotting of the convoluted spectra for direct comparison with inelastic neutron scattering
* full integration into Horace for plotting and comparison with time of flight neutron data, see http://horace.isis.rl.ac.uk
Fitting spin wave spectra
..........................
* possible to fit any parameter of the Hamiltonian
* robust fitting, even when the number of simulated spin wave modes differs from the measured number of modes