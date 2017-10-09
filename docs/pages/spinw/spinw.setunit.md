---
{title: spinw.setunit method, link: spinw.setunit, summary: sets the physical units,
  keywords: sample, sidebar: sw_sidebar, permalink: spinw_setunit, folder: spinw,
  mathjax: 'true'}

---
  
### Syntax
  
`setunit(obj,Name,Value)`
  
### Description
  
`setunit(obj,Name,Value)` sets the physical units of the Hamiltonian.
This includes the magnetic field, exchange interaction, length and
temperature.
  
### Input Arguments
  
`obj`
: [spinw](spinw) object.
  
### Name-Value Pair Arguments
  
`'mode'`
: Type of unit system, defined by one of the following strings:
  * `'AmeVTK'`    Typical units used in neutron/x-ray scattering:
                      [Å, meV, Tesla and Kelvin]
  * `'1'`         No units, all conversion factors are set to 1.
 
### See Also
 
[spinw.unit](spinw_unit)
 

{% include links.html %}
