---
{title: sw_nb, link: sw_nb, summary: returns the bound coherent neutron scattering
    length, keywords: sample, sidebar: sw_sidebar, permalink: sw_nb, folder: swfiles,
  mathjax: 'true'}

---
  
### Syntax
  
`bc = sw_nb(atomname)`
  
### Description
  
`bc = sw_nb(atomname)` returns the bound coherent neutron scattering
length of a given nucleus in fm units. The function reads the stored data
from the `isotope.dat` file.
  
### Input Arguments
  
`atomName`
: String, contains the name of the atom or isotope (e.g. `'13C'` stands
  for the carbon-13 isotope).
  
### Output Arguments
  
`bc`
: Value of the bound coherent neutron scattering length in units of fm.
 

{% include links.html %}
