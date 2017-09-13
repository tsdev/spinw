---
{title: sw_omegasum( ), summary: removes degenerate and ghost magnon modes from spectrum,
  keywords: sample, sidebar: sw_sidebar, permalink: sw_omegasum.html, folder: swfiles,
  mathjax: 'true'}

---
removes degenerate and ghost magnon modes from spectrum
 
spec = [sw_omegasum](sw_omegasum.html)(spec, 'Option1', Value1, ...)
 
It removes the degenerate modes from the dispersion stored in
spectra.omega and sorts omega according to the energy.
 
The degenerate dispersion energies are substituted with NaN values. Be
carefull, after this function [sw_egrid()](sw_egrid.html) won't work properly on spectra.
It doesn't work for spectra with multiple twins.
 
Options:
 
tol       Tolerance, within two energies are considered equal. Default
          value is 1e-5.
zeroint   The minimum intensity value, below the mode is dropped. Default
          value is zero (no modes are dropped due to weak intensity).
emptyval  Value that is assigned to modes, that are removed due to the
          summation. Default value is NaN (good for plotting). Zero can
          be used for further numerical treatmen.
 
See also SPINW.SPINWAVE, SW_EGRID.
 

