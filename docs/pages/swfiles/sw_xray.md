---
{title: sw_xray, link: sw_xray, summary: calculates x-ray scattering cross section,
  keywords: sample, sidebar: sw_sidebar, permalink: sw_xray, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`spectra = sw_xray(spectra,Name,Value)`
  
### Description
  
`spectra = sw_xray(spectra,Name,Value)` calculates x-ray scattering
intensity for non-resonant inelastic x-ray scattering on phonons.
   
  
### Input Arguments
  
`spectra`
: Input structure that contains the displacement-displacement
  correlation function.
  
### Output Arguments
  
`spectra`
: Structure that is same as the input with the following additional
  fields:
  * `param`   Input parameters.
  * `Sperp`   $$S_\perp(i_{mode},\mathbf{Q})$$ x-ray scattering cross
              section, stored in a matrix with dimensions of
              $$[n_{mode n_{hkl}]$$.
  
### See Also
  
[spinw](spinw) \| [spinw.spinwave](spinw_spinwave) \| [sw_neutron](sw_neutron)
 

{% include links.html %}
