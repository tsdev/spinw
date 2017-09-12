---
{title: '@spinw.spinwave( )', summary: calculates dynamical spin-spin correlation
    function using linear spin wave theory, keywords: sample, sidebar: sw_sidebar,
  permalink: spinw_spinwave.html, folder: '@spinw', mathjax: 'true'}

---
calculates dynamical spin-spin correlation function using linear spin wave theory
 
spectra = SPINWAVE(obj, hkl, 'option1', value1 ...)
 
Spin wave dispersion and spin-spin correlation function is calculated at
the reciprocal space points k. The function can deal with arbitrary
magnetic structure and magnetic interactions as well as single ion
anisotropy and magnetic field. Biquadratic exchange interactions are also
implemented, however only for k=0 magnetic structures.
 
If the magnetic ordering wavevector is non-integer, the dispersion is
calculated using a coordinate system rotating from cell to cell. In this
case the spin Hamiltonian has to fulfill this extra rotational symmetry.
 
Some of the code of the function can run faster is mex files are used. To
switch on mex files, use the swpref.setpref('usemex',true) command. For
details see the <a href="matlab:help('sw_mex.m')">sw_mex</a> function.
 
 
Input:
 
obj           Input structure, spinw class object.
hkl           Defines the Q points where the spectra is calculated, in
              reciprocal lattice units, size is [3 nHkl]. Q can be also
              defined by several linear scan in reciprocal space. In this
              case hkl is cell type, where each element of the cell
              defines a point in Q space. Linear scans are assumed
              between consecutive points. Also the number of Q points can
              be specified as a last element, it is 100 by defaults. For
              example: hkl = {[0 0 0] [1 0 0]  50}, defines a scan along
              (h,0,0) from 0 to 1 and 50 Q points are calculated along
              the scan.
 
              For symbolic calculation at a general reciprocal space
              point use sym class input. For example to calculate the
              spectrum along (h,0,0): hkl = [sym('h') 0 0]. To
              do calculation at a specific point do for example
              sym([0 1 0]), to calculate the spectrum at (0,1,0).
 
Options:
 
formfact      If true, the magnetic form factor is included in the
              spin-spin correlation function calculation. The form factor
              coefficients are stored in obj.unit_cell.ff(1,:,atomIndex).
              Default value is false.
formfactfun   Function that calculates the magnetic form factor for given
              Q value. Default value is @sw_mff(), that uses a tabulated
              coefficients for the form factor calculation. For
              anisotropic form factors a user defined function can be
              written that has the following header:
                  F = @formfactfun(atomLabel,Q)
              where the parameters are:
                  F   row vector containing the form factor for every
                      input Q value
                  atomLabel string, label of the selected magnetic atom
                  Q   matrix with dimensions of [3 nQ], where each column
                      contains a Q vector in Angstrom^-1 units.
gtensor       If true, the g-tensor will be included in the spin-spin
              correlation function. Including anisotropic g-tensor or
              different g-tensor for different ions is only possible
              here. Including a simple isotropic g-tensor is possible
              afterwards using the sw_instrument() function.
fitmode       Speedup (for fitting mode only), default is false.
notwin        If true, the spectra of the twins won't be calculated.
              Default is false.
sortMode      The spin wave modes will be sorted if true. Default is
              true.
optmem        Parameter to optimise memory usage. The list of hkl values
              will be cut into optmem number of pieces and will be
              calculated piece by piece to decrease memory usage. Default
              of optmem is zero, when the number of slices are determined
              automatically from the available free memory.
tol           Tolerance of the incommensurability of the magnetic
              ordering wavevector. Deviations from integer values of the
              ordering wavevector smaller than the tolerance are
              considered to be commensurate. Default value is 1e-4.
omega_tol     Tolerance on the energy difference of degenerate modes when
              diagonalising the quadratic form, default is 1e-5.
hermit        Method for matrix diagonalization:
                  true      J.H.P. Colpa, Physica 93A (1978) 327,
                  false     R.M. White, PR 139 (1965) A450.
              Colpa: the grand dynamical matrix is converted into another
                     Hermitian matrix, that will give the real
                     eigenvalues.
              White: the non-Hermitian g*H matrix will be diagonalised,
                     that is not the elegant method.
              Advise:
              Always use Colpa's method, except when small imaginary
              eigenvalues are expected. In this case only White's method
              work. The solution in this case is wrong, however by
              examining the eigenvalues it can give a hint where the
              problem is.
              Default is true.
saveH         If true, the quadratic form of the Hamiltonian is saved. Be
              carefull, it can take up lots of memory. Default is false.
saveV         If true, the matrices that transform the normal magnon
              modes into the magnon modes localized on the spins are
              saved. Be carefull, it can take up lots of memory.
              Default is false.
saveSabp      If true, the dynamical structure factorin the rotating
              frame is saved S'(k,omega). Default is false.
title         Gives a title string to the simulation that is saved in the
              output.
fid           Defines whether to provide text output. Default is defined
              in obj.fid. The possible values:
                  0       No text output is generated.
                  1       Text output in the MATLAB Command Window.
                  fid     File ID provided by the fopen() command, the
                          output is written into the opened file stream.
 
Output:
 
'spectra' is a structure, with the following fields:
omega         Calculated spin wave dispersion, dimensins are
              [nMode nHkl], where nMagExt is the number of magnetic
              atoms in the extended unit cell.
Sab           Dynamical structure factor, dimensins are
              [3 3 nMode nHkl]. Each (:,:,i,j) submatrix contains the
              9 correlation functions: Sxx, Sxy, Sxz, etc. If given,
              magnetic form factor is included. Intensity is in hBar
              units, normalized to the crystallographic unit cell.
H             Quadratic form of the Hamiltonian.
              Only saved if saveH is true.
V             Transformation matrix from the normal magnon modes to the
              magnons localized on spins:
                  x_i = sum_j V_ij * x_j'
              Only saved if saveV is true.
Sabp          Dynamical structure factor in the rotating frame,
              dimensions are [3 3 nMode nHkl], but the number of modes
              are equal to twice the number of magnetic atoms.
formfact      Cell containing the labels of the magnetic ions if form
              factor in included in the spin-spin correlation function.
cmplxBase     The local coordinate system on each magnetic moment is
              defined by the complex magnetic moments:
                  e1 = imag(M/norm(M))
                  e3 = real(M/norm(M))
                  e2 = cross(e3,e1)
 
nMode is the number of magnetic mode. For commensurate structures it is
double the number of magnetic atoms in the magnetic cell/supercell. For
incommensurate structures this number is tripled due to the appearance of
the (Q+/-km) Fourier components in the correlation functions. For every k
points in the following order: (k-km,k,k+km).
 
If several twins exist in the sample, omega and Sab are packaged into a
cell, that contains nTwin number of matrices.
 
hkl           Contains the input Q values, dimensins are [3 nHkl].
hklA          Same Q values, but in reciproc Angstrom units in the
              lab coordinate system, dimensins are [3 nHkl].
incomm        Whether the spectra calculated is incommensurate or not.
obj           The copy of the input obj.
 
Example:
 
tri = sw_model('triAF',1);
sw_plotspec(tri.spinwave({[0 0 0] [1 1 0]}))
 
The above example will calculate and plot the spin wave dispersion of the
triangular lattice antiferromagnet (S=1, J=1) along the [H H 0] direction
in reciprocal space.
 
See also SPINW, SPINW.SPINWAVESYM, SW_MEX, SPINW.POWSPEC.
 
