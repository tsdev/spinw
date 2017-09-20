---
{title: spinw.meanfield method, link: spinw.meanfield, summary: mean field calculation
    of the wave vector dependent susceptibility, keywords: sample, sidebar: sw_sidebar,
  permalink: spinw_meanfield.html, folder: spinw, mathjax: 'true'}

---

### Syntax

`chi = meanfield(obj, hkl,Name,Value)`

### Description



### Input Arguments

`obj`
:nw] object.

`hkl`
:      Defines the Q points where chi is calculated, in reciprocal
       lattice units, size is [3 nHkl]. Q can be also defined by
       several linear scan in reciprocal space. In this case hkl
       is cell type, where each element of the cell defines a
       point in Q space. Linear scans are assumed between
       consecutive points. Also the number of Q points can be
       specified as a last element, it is 100 by defaults. For
       example: hkl = {[0 0 0] [1 0 0]  50}, defines a scan along
       (h,0,0) from 0 to 1 and 50 Q points are calculated along
       the scan.
       For symbolic calculation at a general reciprocal space
       point use sym class input. For example to calculate chi
       along (h,0,0): hkl = [sym('h') 0 0]. To do calculation at a
       specific point do for example sym([0 1 0]), to calculate
       the spectrum at (0,1,0).

### Name-Value Pair Arguments

`'Trel'`
: Relative mean field temperature in Kelvin. An effective
  temperature relative to the mean field critical temperature
  Tc (the most negative eigenvalue of J(Q), the Fourier
  transform of the exchange couplings). Default value is 0,
  which means Tmf = Tc.

`'Tc'`
: Critical temperature, default value is calculated from the
  exchange matrix sampled on the given Q points. If the Q
  points don't contain the point where J(Q) is minimum, the
  automatically determined Tc will be wrong. In this case it
  is recommended to use this option to give the right Tc.

`'formfact'`
: Setting, that determines whether the magnetic form factor
  is included in the spin-spin correlation function
  calculation. Possible values:
      false   No magnetic form factor is applied (default).
      true    Magnetic form factors are applied, based on the
              label string of the magnetic ions, see sw_mff()
              function help.
      cell    Cell type that contains mixed labels and
              numbers for every symmetry inequivalent atom in
              the unit cell, the numbers are taken as
              constants.
  For example: formfact = {0 'MCr3'}, this won't include
  correlations on the first atom and using the form factor of
  Cr3+ ion for the second atom.

`'formfactfun'`
: Function that calculates the magnetic form factor for given
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
          contains a Q vector in Å$$^{-1}$$ units.

### Output Arguments

chi           Structure with the following fields:
Tc            Critical temperature in Kelvin, determined from the minimum
              eigenvalue of J(Q) sampled on the Q points given in the
              input.
Tmf           Mean field temperature, Tmf = Trel + Tc.
...

