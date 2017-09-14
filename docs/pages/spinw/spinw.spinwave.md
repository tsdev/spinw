---
{title: spinw.spinwave method, link: spinw.spinwave, summary: calculates spin correlation
    function using linear spin wave theory, keywords: sample, sidebar: sw_sidebar,
  permalink: spinw_spinwave.html, folder: spinw, mathjax: 'true'}

---
 
* * *
`spectra = spinwave(obj, hkl, 'option1', value1 ...)`
* * *
 
Spin wave dispersion and spin-spin correlation function is calculated at
the reciprocal space points $$k$$. The function can deal with arbitrary
magnetic structure and magnetic interactions as well as single ion
anisotropy and magnetic field. Biquadratic exchange interactions are also
implemented, however only for $$k=0$$ magnetic structures.
 
If the magnetic ordering wavevector is non-integer, the dispersion is
calculated using a coordinate system rotating from cell to cell. In this
case the spin Hamiltonian has to fulfill this extra rotational symmetry
which is not checked programatically.
 
Some of the code of the function can run faster if mex files are used. To
switch on mex files, use the `swpref.setpref('usemex',true)` command. For
details see the [sw_mex](sw_mex.html) and [swpref.setpref](swpref_setpref.html) functions.
 
 
### Input
 
`obj`
: [spinw](spinw.html) object.
  
`hkl`
: Defines the $$Q$$ points where the spectra is calculated, in reciprocal
  lattice units, size is $$[3\times n_{hkl}]$$. $$Q$$ can be also defined by
  several linear scan in reciprocal space. In this case `hkl` is cell type,
  where each element of the cell defines a point in $$Q$$ space. Linear scans
  are assumed between consecutive points. Also the number of $$Q$$ points can
  be specified as a last element, it is 100 by defaults. 
    
  For example to define a scan along $$(h,0,0)$$ from $$h=0$$ to $$h=1$$ using
  200 $$Q$$ points the following input should be used:
  ```matlab
  hkl = {[0 0 0] [1 0 0]  50}
  ```
 
  For symbolic calculation at a general reciprocal space point use `sym`
  type input. 
 
  For example to calculate the spectrum along $$(h,0,0)$$ use:
  ```matlab
  hkl = [sym('h') 0 0]
  ```
  To calculate spectrum at a specific $$Q$$ point symbolically, e.g. at
  $$(0,1,0)$$ use:
  ```matlab
  hkl = sym([0 1 0])
  ```
 
### Options
 
`formfact`
: If true, the magnetic form factor is included in the spin-spin
  correlation function calculation. The form factor coefficients are
  stored in `obj.unit_cell.ff(1,:,atomIndex)`. Default value is `false`.
 
`formfactfun`
: Function that calculates the magnetic form factor for given $$Q$$ value.
  value. Default value is `@sw_mff`, that uses a tabulated coefficients
  for the form factor calculation. For anisotropic form factors a user
  defined function can be written that has the following header:
  ```matlab
  F = formfactfun(atomLabel,Q)
  ```
  where the parameters are:
  * `F`           row vector containing the form factor for every input 
                  $$Q$$ value
  * `atomLabel`   string, label of the selected magnetic atom
  * `Q`           matrix with dimensions of $$[3\times n_Q]$$, where each
                  column contains a $$Q$$ vector in $$Å^{-1}$$ units.
 
`gtensor`
: If true, the g-tensor will be included in the spin-spin correlation
  function. Including anisotropic g-tensor or different
  g-tensor for different ions is only possible here. Including a simple
  isotropic g-tensor is possible afterwards using the [sw_instrument](sw_instrument.html)
  function.
 
`fitmode`
: If `true`, function is optimized for multiple consecutive calls (e.g. 
  the output spectrum won't contain the copy of `obj`), default is
  `false`.
 
`notwin`
: If `true`, the spectra of the twins won't be calculated. Default is
`false`.
 
`sortMode`
: If `true`, the spin wave modes will be sorted. Default is `true`.
 
`optmem`
: Parameter to optimise memory usage. The list of hkl values will be cut
  into `optmem` number of pieces and will be calculated piece by piece to
  decrease peak memory usage. Default value is 0, when the number
  of slices are determined automatically from the available free memory.
 
`tol`
: Tolerance of the incommensurability of the magnetic ordering wavevector.
  Deviations from integer values of the ordering wavevector smaller than
  the tolerance are considered to be commensurate. Default value is
  $$10^{-4}$$.
 
`omega_tol`
: Tolerance on the energy difference of degenerate modes when
  diagonalising the quadratic form, default value is $$10^{-5}$$.
 
`hermit`
: Method for matrix diagonalization with the following logical values:
  
  * `true`    using Colpa's method (for details see [J.H.P. Colpa, Physica 93A (1978) 327](http://www.sciencedirect.com/science/article/pii/0378437178901607)),
              the dynamical matrix is converted into another Hermitian
              matrix, that will give the real eigenvalues.
  * `false`   using the standard method (for details see [R.M. White, PR 139 (1965) A450](https://journals.aps.org/pr/abstract/10.1103/PhysRev.139.A450))
              the non-Hermitian $$\mathcal{g}\times \mathcal{H}$$ matrix
              will be diagonalised, which is computationally less
              efficient. Default value is `true`.
 
{% include note.html content="
  Always use Colpa's method, except when imaginary eigenvalues are
  expected. In this case only White's method work. The solution in this
  case is wrong, however by examining the eigenvalues it can give a hint
  where the problem is." %}
                
`saveH`
: If true, the quadratic form of the Hamiltonian is also saved in the
  output. Be carefull, it can take up lots of memory. Default value is
  `false`.
 
`saveV`
: If true, the matrices that transform the normal magnon modes into the
  magnon modes localized on the spins are also saved into the output. Be
  carefull, it can take up lots of memory. Default value is `false`.
 
`saveSabp`
: If true, the dynamical structure factor in the rotating frame
  $$S'(k,\omega)$$ is saved. Default value is `false`.
 
`title`
: Gives a title string to the simulation that is saved in the output.
 
`fid`
: Defines whether to provide text output. Default value is defined in
  `obj.fid`. The possible values are:
                    
  * `0`   No text output is generated.
  * `1`   Text output in the MATLAB Command Window.
  * `fid` File ID provided by the `fopen` command, the output is written
          into the opened file stream.
 
### Output
 
`spectra`
: structure, with the following fields:
  
  * `omega`   Calculated spin wave dispersion with dimensions of
              $$[n_{mode}\times n_{hkl}]$$.
  * `Sab`     Dynamical structure factor with dimensins of
              $$[3\times 3\times n_{mode}\times n_{hkl}]$$. Each
              `(:,:,i,j)` submatrix contains the 9 correlation functions
              $$S^{xx}$$, $$S^{xy}$$, $$S^{xz}$$, etc. If given, magnetic form
              factor is included. Intensity is in ħ units, normalized
              to the crystallographic unit cell.
  * `H`       Quadratic form of the Hamiltonian. Only saved if `saveH` is
              true.
  * `V`       Transformation matrix from the normal magnon modes to the
              magnons localized on spins using the following:
              $$x_i = \sum_j V_{ij} \times x_j'$$
              Only saved if `saveV` is true.
  * `Sabp`    Dynamical structure factor in the rotating frame,
              dimensions are $$[3\times 3\times n_{mode}\times n_{hkl}]$$,
              but the number of modes are equal to twice the number of
              magnetic atoms.
  * `formfact`  Cell containing the labels of the magnetic ions if form
              factor in included in the spin-spin correlation function.
  * `cmplxBase` The local coordinate system on each magnetic moment is
              defined by the complex magnetic moments:
              $$\begin{align}  e_1 &= \Im(\hat{M})\\
                              e_3 &= Re(\hat{M})\\
                              e_2 &= e_3\times e_1
              \end{align}$$
 
  * `hkl`     Contains the input $$Q$$ values, dimensions are $$[3\times n_{hkl}]$$.
  * `hklA`    Same $$Q$$ values, but in $$Å^{-1}$$ unit, in the
              lab coordinate system, dimensins are $$[3\times n_{hkl}]$$.
  * `incomm`  Logical value, tells whether the calculated spectra is
              incommensurate or not.
  * `obj`     The copy (clone) of the input `obj`, see [spinw.copy](spinw_copy.html).
 
The number of magnetic modes (labelled by `nMode`) for commensurate
structures is double the number of magnetic atoms in the magnetic cell.
For incommensurate structures this number is tripled due to the
appearance of the $$(Q\pm k_m)$$ Fourier components in the correlation
functions. For every $$Q$$ points in the following order:
$$(Q-k_m,Q,Q+k_m)$$.
 
If several twins exist in the sample, `omega` and `Sab` are packaged into
a cell, that contains $$n_{twin}$$ number of matrices.
 
### Example
 
```matlab
  tri = sw_model('triAF',1);
  sw_plotspec(tri.spinwave({[0 0 0] [1 1 0]}))
```
 
The above example will calculate and plot the spin wave dispersion of the
triangular lattice antiferromagnet ($$S=1$$, $$J=1$$) along the $$(h,h,0)$$
direction in reciprocal space.
 
### See also
[spinw](spinw.html), [spinw.spinwavesym](spinw_spinwavesym.html), [sw_mex](sw_mex.html) and [spinw.powspec](spinw_powspec.html)
 

