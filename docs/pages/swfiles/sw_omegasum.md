---
{title: sw_omegasum( ), link: sw_omegasum, summary: removes degenerate and ghost magnon
    modes from spectrum, keywords: sample, sidebar: sw_sidebar, permalink: sw_omegasum.html,
  folder: swfiles, mathjax: 'true'}

---

### Syntax

` `

### Description

spectra.omega and sorts omega according to the energy.
 
The degenerate dispersion energies are substituted with NaN values. Be
carefull, after this function sw_egrid() won't work properly on spectra.
It doesn't work for spectra with multiple twins.
 

### Name-Value Pair Arguments

% `tol`
:   Tolerance, within two energies are considered equal. Default

% `value`
:s 1e-5.

% `zeroint`
:   The minimum intensity value, below the mode is dropped. Default

% `value`
:s zero (no modes are dropped due to weak intensity).

% `emptyval`
:l  Value that is assigned to modes, that are removed due to the

% `summation.`
:on. Default value is NaN (good for plotting). Zero can

% `be`
: for further numerical treatmen.

### See Also

[spinw.spinwave](spinw_spinwave.html) and [sw_egrid](sw_egrid.html)

