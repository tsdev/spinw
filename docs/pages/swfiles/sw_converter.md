---
{title: sw_converter, link: sw_converter, summary: converts energy and momentum units
    for a given particle, keywords: sample, sidebar: sw_sidebar, permalink: sw_converter.html,
  folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`out = sw_converter(value, unitin, unitout, {particlename}) `
  
### Description
  
`out = sw_converter(value, unitin, unitout, {particlename}) `
  
### Input Arguments
  
`particleName`
: Name of the particle:
      'neutron'   default
      'proton'
      'electron'
      'photon'
  
`value`
: Numerical input value, can be arbitrary matrix.
  
`unitIn`
: Units of the input value:
      'A-1'       momentum in Å$$^{-1}$$.
      'k'         -||-
      'Å'  wavelength in Å.
      'lambda'    -||-
      'A'         -||-
      'Kelvin'    temperature in Kelvin.
      'K'         -||-
      'mps'       speed in m/s.
      'J'         energy in Joule.
      'meV'       energy in meV.
      'THz'       frequency in Thz.
      'cm-1'      2*pi/lambda in cm$$^{-1}$$.
      'eV'        energy in eV.
      'fs'        wave period time in fs.
      'ps'        wave period time in ps.
      'nm'        wavelength in nm.
      'um'        wavelength in um.
  
`unitOut`
: Units of the output value, same options as for unitIn.
 

{% include links.html %}
