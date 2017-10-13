---
{title: sw_converter, link: sw_converter, summary: converts energy and momentum units
    for a given particle, keywords: sample, sidebar: sw_sidebar, permalink: sw_converter,
  folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`out = sw_converter(value, unitIn, unitOut)`
  
`out = sw_converter(value, unitIn, unitOut, particleName)`
 
### Description
  
`out = sw_converter(value, unitin, unitout)` will convert momentum and
energy values assuming neutron as a particle.
  
`out = sw_converter(value, unitin, unitout,particleName)` will convert
momentum and energy values for a given particle, such as neutron, photon,
etc.
 
### Example
 
Calculate the energy of a neutron (in meV) which has a wavelength of
5 Å:
 
```matlab
sw_converter(5,'A','meV')
```
*Output*
```
    3.2722
```
 
 
Calculate the wavelength of X-ray in Å that has 7.5 keV energy:
 
```matlab
sw_converter(7.5,'keV','A','photon')
```
*Output*
```
    1.6531
```
 
 
### Input Arguments
  
`value`
: Numerical input value, can be scalar or matrix with arbitrary
  dimensions.
 
`unitIn`
: Units of the input value, one of the following string:
  * `'A-1'`        momentum in Å$$^{-1}$$,
  * `'k'`          momentum in Å$$^{-1}$$,
  * `'Angstrom'`   wavelength in Å,
  * `'lambda'`     wavelength in Å,
  * `'A'`          wavelength in Å,
  * `'K'`          temperature in Kelvin,
  * `'m/s'`        speed in m/s,
  * `'J'`          energy in Joule,
  * `'meV'`        energy in meV,
  * `'eV'`         energy in eV,
  * `'keV'`        energy in keV,
  * `'THz'`        frequency in Thz,
  * `'cm-1'`       $$2\pi/\lambda$$ in cm$$^{-1}$$,
  * `'fs'`         wave period time in fs,
  * `'ps'`         wave period time in ps,
  * `'nm'`         wavelength in nm,
  * `'um'`         wavelength in $$\mu$$m.
  
`unitOut`
: Units of the output value, same strings are accepted as for `unitIn`.
 
`particleName`
: String, the name of the particle, one of the following values: 
  `'neutron'` (default), `'proton'`, `'electron'`, `'photon'`, `'xray'`,
  `'light'`.
 

{% include links.html %}
