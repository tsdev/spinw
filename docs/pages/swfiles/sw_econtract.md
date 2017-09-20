---
{title: sw_econtract( ), link: sw_econtract, summary: 'converts (Q,E) values to Q
    values for diffraction instrument', keywords: sample, sidebar: sw_sidebar, permalink: sw_econtract.html,
  folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`Qm = sw_econtract(Q,Name,Value)`
  
### Description
  
`Qm = sw_econtract(Q,Name,Value)` converts $$(Q,E)$$ values in the phase
space into $$Q$$ values as it would appear when measuring it with via
neutron diffraction.
  
### Input Arguments
  
`Q`
: Input values in reciprocal space in the scattering plane in
  Å$$^{-1}$$ units, dimensions are $$[2\times n_Q]$$.
  
### Name-Value Pair Arguments
  
`omega`
: Energy transfer value in meV, default value is zero.
  
`lambda`
: Wavelength of the incident neutron beam in Å.
  
`ki`
: Momentum of the incidend neutron beam in Å$$^{-1}$$, alternative
  input to `lambda`.
  
`sense`
: Scattering sense:
 
  * `1`  detectors are on the right hand side from the incident beam direction, default.
  * `-1` detectors are on the left hand side from the incident beam direction.
  
### See Also
  
[sw_converter](sw_converter.html)
 

