---
{title: spinw.spinwavesym method, link: spinw.spinwavesym, summary: calculates symbolic
    spin wave dispersion, keywords: sample, sidebar: sw_sidebar, permalink: spinw_spinwavesym.html,
  folder: spinw, mathjax: 'true'}

---
 
spectra = SPINWAVESYM(obj, 'option1', value1 ...)
 
Symbolic spin wave dispersion is calculated as a function of reciprocal
space points. The function can deal with arbitrary magnetic structure and
magnetic interactions as well as single ion anisotropy and magnetic
field. Biquadratic exchange interactions are also implemented, however
only for k=0 magnetic structures.
 
If the magnetic ordering wavevector is non-integer, the dispersion is
calculated using a coordinate system rotating from cell to cell. In this
case the Hamiltonian has to fulfill this extra rotational symmetry.
 
The method works for incommensurate structures, however the calculated
omega dispersion does not contain the omega(k+/-km) terms that has to be
added manually.
 
The method for matrix diagonalization is according to R.M. White, PR 139
(1965) A450. The non-Hermitian g*H matrix will be diagonalised.
 
Input:
 
obj           Input structure, spinw class object.
 
Options:
 
hkl           Symbolic definition of q vector. Default is the general Q
              point:
                  hkl = [sym('h') sym('k') sym('l')]
eig           If true the symbolic Hamiltonian is diagonalised. For large
              matrices (many magnetic atom per unit cell) this might be
              impossible. Set 'eig' to false to output only the quadratic
              Hamiltonian. Default is true.
vect          If true the eigenvectors are calculated. Default is false.
tol           Tolerance of the incommensurability of the magnetic
              ordering wavevector. Deviations from integer values of the
              ordering wavevector smaller than the tolerance are
              considered to be commensurate. Default value is 1e-4.
norm          Whether to produce the normalized eigenvectors. It can be
              impossible for large matrices, in that case set it to
              false. Default is true.
title         Gives a title string to the simulation that is saved in the
              output.

