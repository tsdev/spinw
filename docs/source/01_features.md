---
title: Features of SpinW
keywords: homepage
sidebar: sw_sidebar
permalink: features.html
summary: Quick list of features
folder: documentation
mathjax: true
---

SpinW is a Matlab library that can define spin Hamiltonians including general quadratic exchange interactions, single ion anisotropy and external magnetic field, applying space group symmetries and solve the resulting equation of motion.

### Summary

In short SpinW can solve the following spin Hamiltonian using classical and quasi classical numerical methods:

$\mathcal{H} = \sum_{i,j} \textbf{S}_i J_{ij}\textbf{S}_j + \sum_{i,j}\textbf{S}_i A_i\textbf{S}_j + \textbf{B}\sum_i g_i\textbf{S}_i$

where $\textbf{S}_i$ are spin vector operators, $J_{ij}$ are $[3\times 3]$ matrices describing pair coupling between spins, $A_{ij}$ are $[3\times 3]$ anisotropy matrices, $\textbf{B}$ is external magnetic field and $g_i$ is the g-tensor.

### Crystal structures

* definition of crystal structure with arbitrary unit cell, using space group or symmetry operators
* definition of non-magnetic atoms and magnetic atoms with arbitrary moment size
* publication quality 3D plotting of crystal structures (atoms, labels, axes, surrounding polyhedron, anisotropy ellipsoids, DM vector, etc.)

### Magnetic structures

* definition of magnetic structures using complex magnetization vectors
* representation of incommensurate structures using rotating frame coordinate system
* generation of magnetic structures on a magnetic supercell
* 3D visualization of magnetic structures

### Magnetic interactions

* simple assignment of bonds based on length
* arbitrary quadratic exchange interactions are allowed (including DM, etc.)
* arbitrary single ion anisotropy tensor (easy-plane, easy-axis, etc.)
* Zeeman energy in homogeneous magnetic field including arbitrary g-tensor
* calculation of symmetry allowed elements of the above tensors based on the crystallographic space group

### Optimization of magnetic structures

* minimization of the classical energy assuming single-k magnetic structure for fast and simple solution for ground state magnetic structure
* simulated annealing using the Metropolis algorithm on an arbitrary large magnetic supercell
* calculating properties in thermodynamical equilibrium (heat capacity, magnetic susceptibility, etc.)
* calculation of the magnetic structure factor
* simulation of magnetic neutron diffraction and diffuse scattering

### Linear spin wave theory

* simulation of magnetic excitations in general commensurate and incommensurate magnetic structures using linear spin-wave theory
* calculation of spin wave dispersion, spin-spin correlation functions
* calculation of neutron scattering cross section for unpolarized neutrons including the magnetic form factor
* calculation of polarized neutron scattering cross sections
* possible to include different moment sizes for different magnetic atoms
* calculation of powder averaged spin wave spectrum

### Plotting of spin wave spectrum

* plotting of dispersions and correlation functions
* calculation and plotting of the convoluted spectra for direct comparison with inelastic neutron scattering
* full integration into Horace for plotting and comparison with time of flight neutron data, see http://horace.isis.rl.ac.uk

### Fitting spin wave spectra

* possible to fit any parameter in the Hamiltonian
* robust fitting, even when the number of simulated spin wave modes differs from the measured number of modes

