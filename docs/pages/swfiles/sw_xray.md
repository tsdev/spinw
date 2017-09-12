---
{title: sw_xray( ), summary: calculates X-ray scattering intensity for phonon spectrum,
  keywords: sample, sidebar: sw_sidebar, permalink: swfiles_sw_xray.html, folder: swfiles,
  mathjax: 'true'}

---
calculates X-ray scattering intensity for phonon spectrum
 
spectra = SW_XRAY(spectra, 'Option1', Value1' ...)
 
It calculates X-ray scattering intensity for inelastic X-ray scattering
on phonons.
 
Input:
 
spectra   Input structure that contains the displacement-displacement
          correlation function.
 
Output:
 
the spectra output has the following additional fields:
param     Input parameters.
 
Sperp     Sperp(mode,Q) X-ray scattering cross section, dimensions are
          [nMode nHkl].
formfact      Setting, that determines whether the magnetic form factor
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
 
If several domains exist in the sample, Sperp is packaged into a cell,
that contains nTwin number of matrices.
 
See also SPINW, SPINW.SPINWAVE, SW_NEUTRON.
 
