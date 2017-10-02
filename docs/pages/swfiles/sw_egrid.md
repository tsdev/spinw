---
{title: sw_egrid, link: sw_egrid, summary: calculates energy bins of a spectrum, keywords: sample,
  sidebar: sw_sidebar, permalink: sw_egrid, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`spectra = sw_egrid(spectra,Name,Value)`
  
### Description
  
`spectra = sw_egrid(spectra,Name,Value)` takes a calculated spectrum that
contains $$S^{\alpha\beta}(Q,\omega)$$ and converts it into an intensity
map `I(i,j)` via binning the energy values and selecting a given
component of the $$9\times 9$$ spin-spin correlation matrix. For example by
default (setting the `component` parameter to `'Sperp'`) it selects the
neutron scattering cross section via calculating the following quantity:
 
  $$S_\perp(Q,\omega)=\sum_{\alpha\beta}(1-\hat{q}^\alpha\hat{q}^\beta)\cdot S^{\alpha\beta}(Q,\omega)$$
   
  
### Examples
  
The line will create an energy bin, with steps of 0.1 and bins the
spin-spin correlation function. Two different matrices will be
calculated, first using the sum of the $$S^{xx}$$ and $$S^{yy}$$ components, second
will contain the $$S^{zz}$$ component of the correlation function.
 
```matlab
tri = sw_model('triAF',1)
spectra = tri.spinwave({[0 0 0] [1 1 0] 501})
E = linspace(0,5,501)
spectra = sw_egrid(spectra,'component',{'Sxx+Syy' 'Szz'},'Evect',E)
figure
sw_plotspec(spectra,'mode','color','axLim',[0 0.5],'dE',0.2)
```
 
{% include image.html file="generated/sw__1.png" alt="sw_plotspec(spectra,'mode','color','axLim',[0 0.5],'dE',0.2)" %}
 
### Input Arguments
  
`spectra`
: Input structure, contains spin-spin correlation functions. Supported
  inputs are produced by [spinw.spinwave](spinw_spinwave), [spinw.powspec](spinw_powspec) and
  [spinw.scga](spinw_scga).
  
### Name-Value Pair Arguments
  
`'component'`
: A string that Selects a correlation function component that will be
  binned. The possible values are:
  * `'Sperp'` bins the magnetic neutron scattering intensity
    (the $$\langle S_\perp S_\perp\rangle$$ expectation value). Default.
  * `'Sab'`   bins the selected components of the spin-spin
              correlation function. Letter `a` and `b` can be `x`,
              `y` or `z`. For example: `'Sxx'` will convolute the
              $$S^{xx}(Q,\omega)$$ component of the correlation function with the
              dispersion. Here the $$xyz$$ is the standard coordinate system.
  *`'Mab'`    bins the selected components of the spin-spin
              correlation function in the Blume-Maleev coordinate system.
              Letter `a` and `b` can be `x`, `y` or `z`. For example:
              `'Mxx'` will convolute the `xx` component of the
              correlation function with the dispersion.
  * `'Pab'`   bins the selected component of the polarisation
              matrix. Letter `a` and `b` can be `x`, `y` or `z`. For
              example: `'Pyy'` will convolute the `yy` component of
              the polarisation matrix with the dispersion. The
              coordinates used are in the Blume-Maleev coordinate
              system, see below.
  * `'Pa'`    bins the intensity of the calculated polarised
              neutron scattering, with inciden polarisation of
              `Pa` where letter `a` can be `x`, `y` or `z`. For example:
              `'Py'` will convolute the scattering intensity
              simulated for incident polarisation $$P_i\|y$$. The
              used coordinates are in the Blume-Maleev coordinate
              system.
  * `'fName'` where `fName` is one of the field names of the input
              structure spectra. This field should contain a
              matrix with dimensions of $$[n_{mode}\times n_{hkl}]$$.
 
  Any linear combination of the above are allowed, for example:
  `'Sxx+2*Syy'` will bin the linear combination of the `xx` component of
  the spin-spin correlation function with the `yy` component.
  Several cross section can be convoluted and stored
  independently, if component is a cell array containing strings
  each containing any linear combination of cross sections as
  above, the cell array needs to have size $$[1\times n_{cell}]$$, for
  example `{'Sxx' 'Syy' 'Szz'}`.
  
`'Evect'`
: Row vector that defines the center/edge of the energy bins of the
  calculated output, number of elements is $$n_E$$. The energy units
  are defined by the [spinw.unit](spinw_unit) property. Default
  value is an edge bin: `linspace(0,1.1*maxOmega,501)`.
  
`'binType'`
: String, determines the type of bin give, possible options:
  * `'cbin'`      Center bin, the center of each energy bin is given.
  * `'ebin'`      Edge bin, the edges of each bin is given.
  Default value is `'ebin'`.
  
`'T'`
: Temperature, used to calculate the Bose factor in units
  depending on the Boltzmann constant stored in [spinw.unit](spinw_unit). Default
  temperature is taken from `obj.single_ion.T`. The Bose factor is
  included in `swConv` field of the output.
  
`'sumtwin'`
: If true, the spectra of the different twins will be summed
  together weighted with the normalized volume fractions, see
  [spinw.twin](spinw_twin). Default value is true.
  
`'modeIdx'`
: Select certain spin wave modes from the $$2*n_{magatom}$$ number of
  modes to include in the output. Default value is `1:2*nMagAtom` to
  include all modes.
  
`'epsilon'`
: Error limit, used to determine whether a given energy bin is
  uniform or not. Default value is $$10^{-5}$$.
  
`'autoEmin'`
: Due to the finite numerical precision, the spin wave energies
  can contain small imaginary values. These can ruin the
  convoluted spectrum at low energies. To improve the spectrum,
  the lowest energy bin should start above the imaginary part of
  the spin wave energy. If `'autoEmin'` is set to `true`, it
  calculates the bottom of the first energy bin automatically and
  overwrites the given value. Only works if the input energy bin
  starts with zero. Default value is `false`.
  
`'imagChk'`
: Checks whether the imaginary part of the spin wave dispersion is
  smaller than the energy bin size. Default value is true.
  
{% include note.html content=" The Blume-Maleev coordinate system is a cartesian coordinate
system with $$x_{BM}$$, $$y_{BM}$$ and $$z_{BM}$$ basis vectors defined as:
<br> $$x_{BM}$$    parallel to the momentum transfer $$Q$$,
<br> $$y_{BM}$$    perpendicular to $$x_{BM}$$ in the scattering plane,
<br> $$z_{BM}$$    perpendicular to the scattering plane.
" %}
  
### Output Arguments
  
`spectra` same as the input `spectra` plus additions fields:
  
`swConv`
: Stores the selected cross section binned in energy in a matrix with
  dimensions of $$[n_E\times n_{hkl}]$$. Includes the Bose factor.
  
`swInt`
: Stores the selected cross sections for every mode in a matrix with
  dimensions of $$[n_{mode}\times n_{hkl}]$$.
  
`T`
: Input temperature.
  
`component`
: Cell that contains the input component selector strings.
  
`Evect`
: Input energy bin vector, defines the energy bin **edge** positions
  (converted from the given bin centers if necessary).
  
`param`
: All the input parameters.
 
If `'component'` parameter is a cell array or the spectra of multiple
twins are convoluted separately, swConv and swInt will be a cell that
packages the matrices corresponding to each component/twin. The
dimensions of the cell are $$[n_{conv}\times n_{twin}]$$.
  
### See Also
  
[spinw.spinwave](spinw_spinwave) \| [sw_neutron](sw_neutron)
 

{% include links.html %}
