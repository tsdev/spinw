---
{title: setunit( ), keywords: sample, summary: sets the physical units, sidebar: sw_sidebar,
  permalink: '@spinw_setunit.html', folder: '@spinw', mathjax: 'true'}

---
  sets the physical units
 
  SETUNIT(obj, 'option1', value1 ...)
 
  Input:
 
  obj       SpinW object.
 
  Options:
 
  mode      Type of unit system, defined by one of the following strings:
                'AmeVTK'    Typical units used in neutron/xray scattering:
                                [Angstrom, meV, Tesla and Kelvin]
                '1'         No units, all conversion factors are set to 1.
 
