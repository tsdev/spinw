---
{title: sw_instrument, link: sw_instrument, summary: convolutes spectrum with different
    functions, keywords: sample, sidebar: sw_sidebar, permalink: sw_instrument.html,
  folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`spectra = sw_instrument(spectra,Name,Value)`
  
### Description
  
`spectra = sw_instrument(spectra,Name,Value)` can convolute an energy
binned spectrum with different energy resolution functions and add other
effects that introduced by measurement (such as the kinematic limit for
neutron scattering, finite momentum resolution or finite detector
coverage).
   
  
### Name-Value Pair Arguments
  
`'dE'`
: Convolutes the spectrum with a Gaussian in energy, where the width is
  defined by the FWHM value. The accepted values are:
  * *string*   File name, that contains the FWHM energy
               resolution values as a function of energy
               transfer. The file has to contain two columns,
               first is the energy values, the second is the
               FWHM resolution at the given energy transfer
               value, see [sw_res](sw_res.html) function for details.
  * *number*   Constant FWHM energy resolution given by the number.
  * *matrix*   Dimensions of $$[N\times 2]$$, where the first column contains the
               energy transfer values, second column contains
               the FWHM resolution values. These discrete values will
               be fitted using a polynomial with a fixed
               degree, see [sw_res](sw_res.html) for details.
  * *function* Function handle of a resolution function
  with the following header:
  ```matlab
  E_FWHM = res_fun(E)
  ```
  where `E_FWHM` is the FWHM energy resolution and `E` is the energy transfer
  value.
  
`'func'`
: Shape of the energy resolution function if different from Gaussian.
  For details see [sw_resconv](sw_resconv.html).
  
`'polDeg'`
: Degree of the polynomial that is fitted to the discrete energy 
  resolution data. Only used if `dE` is a matrix of string. Default value
  is 5.
  
`'dQ'`
: Momentum transfer resolution of the instrument, FWHM is
  given in Å$$^{-1}$$ units by default, unless different units
  are defined in [spinw.unit](spinw_unit.html). Default value is 0 for no convolution.
  
`'thetaMin'`
: Minimum scattering angle in °, default value is 0. Can be only
  applied if one of the `ki`, `Ei`, `kf` or `Ef` parameters is defined.
  
`'plot'`
: If the resolution is read from file and plot option is
  true, the energy dependent resolution values together with the
  polynomial fit will be plotted in a new figure. Default value is
  `true`.
  
`'norm'`
: If true, the data is normalized to mbarn units. Default is
  false. If no g-tensor is included in the spin wave
  calculation, $$g = 2$$ will be assumed for the conversion.
  
`'useRaw'`
: If `false`, the already modified `spectra.swConv` field is
  modified further instead of the original powder spectrum
  stored in `spectra.swRaw`. Default value is `true`.
  
For simulating the effect of the neutron kinematic limit or the finite 
detector coverage of a neutron spectrometer one of the following
parameter has to be given. The unit of these quantities is defined in
[spinw.unit](spinw_unit.html) with default momentum unit of Å$$^{-1}$$ and energy
unit of meV.
 
`'ki'`
: Fixed momentum of the incident neutrons.
  
`'Ei'`
: Fixed energy of the incident neutrons.
  
`'kf'`
: Fixed final momentum of the neutrons.
  
`'Ef'`
: Fixed final energy of the neutrons.
  
  
### Output Arguments
  
`spectra`
: Struct variable, same as input with following additional fields:
* `norm`      `true`, if the spectrum is normalised to mbarn units.
* `ki`        Fixed incident neutron wave vector if defined in the input.
* `kf`        Fixed final neutron wave vector if defined in the input.
* `dE`        Energy resolution polynomial as given in the input.
* `dQ`        FWHM of the momentum resolution.
* `swRaw`     Original simulated data before the application of
              `sw_instrument`.
  
### See Also
  
[polyfit] \| [polyval] \| [sw_res](sw_res.html) \| [sw_resconv](sw_resconv.html)
 
[FWHM]: Full Width at Half Maximum
 

{% include links.html %}
